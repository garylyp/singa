/************************************************************
*
* Licensed to the Apache Software Foundation (ASF) under one
* or more contributor license agreements.  See the NOTICE file
* distributed with this work for additional information
* regarding copyright ownership.  The ASF licenses this file
* to you under the Apache License, Version 2.0 (the
* "License"); you may not use this file except in compliance
* with the License.  You may obtain a copy of the License at
*
*   http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing,
* software distributed under the License is distributed on an
* "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
* KIND, either express or implied.  See the License for the
* specific language governing permissions and limitations
* under the License.
*
*************************************************************/

#include "./cifar10.h"
#include "singa/model/feed_forward_net.h"
#include "singa/model/optimizer.h"
#include "singa/model/metric.h"
#include "singa/utils/channel.h"
#include "singa/utils/string.h"
#include "singa/io/snapshot.h"
#ifdef MY_FILE_READER
#include "file_reader.h"
#endif
namespace singa
{
// currently supports 'cudnn' and 'singacpp'
#ifdef USE_CUDNN
	const std::string engine = "cudnn";
#else
	const std::string engine = "singacpp";
#endif // USE_CUDNN

	LayerConf GenConvConf(string name, int nb_filter, int kernel, int stride, int pad)
	{
		LayerConf conf;
		conf.set_name(name);
		conf.set_type(engine + "_convolution");
		ConvolutionConf *conv = conf.mutable_convolution_conf();
		conv->set_num_output(nb_filter);
		conv->add_kernel_size(kernel);
		conv->add_stride(stride);
		conv->add_pad(pad);
		conv->set_bias_term(true);

		ParamSpec *wspec = conf.add_param();
		wspec->set_name(name + "_weight");
		auto wfill = wspec->mutable_filler();
		wfill->set_type("Gaussian");
		wfill->set_std(sqrt(2.0f/(nb_filter*9.0f)));

		ParamSpec *bspec = conf.add_param();
		bspec->set_name(name + "_bias");
		auto bfill = bspec->mutable_filler();
  		bfill->set_value(0.0f);
		return conf;
	}

	LayerConf GenPoolingConf(string name, bool max_pool, int kernel, int stride, int pad)
	{
		LayerConf conf;
		conf.set_name(name);
		conf.set_type(engine + "_pooling");
		PoolingConf *pool = conf.mutable_pooling_conf();
		pool->set_kernel_size(kernel);
		pool->set_stride(stride);
		pool->set_pad(pad);
		if (!max_pool)
			pool->set_pool(PoolingConf_PoolMethod_AVE);
		return conf;
	}

	LayerConf GenReLUConf(string name)
	{
		LayerConf conf;
		conf.set_name(name);
		conf.set_type(engine + "_relu");
		return conf;
	}

	LayerConf GenDenseConf(string name, int num_output)
	{
		LayerConf conf;
		conf.set_name(name);
		conf.set_type("singa_dense");
		DenseConf *dense = conf.mutable_dense_conf();
		dense->set_num_output(num_output);

		ParamSpec *wspec = conf.add_param();
		wspec->set_name(name + "_weight");
		wspec->set_lr_mult(1);
		wspec->set_decay_mult(0);
		auto wfill = wspec->mutable_filler();
		wfill->set_type("Gaussian");
		wfill->set_std(0.0005);

		ParamSpec *bspec = conf.add_param();
		bspec->set_name(name + "_bias");
		bspec->set_lr_mult(2);
		bspec->set_decay_mult(0);

		return conf;
	}

	LayerConf GenSoftmax(string name)
	{
		LayerConf conf;
		conf.set_name(name);
		conf.set_type(engine + "_softmax");
		SoftmaxConf *softmax = conf.mutable_softmax_conf();
		return conf;
	}

	LayerConf GenFlattenConf(string name)
	{
		LayerConf conf;
		conf.set_name(name);
		conf.set_type("singa_flatten");
		return conf;
	}

	LayerConf GenDropoutConf(string name, float dropout_ratio)
	{
		LayerConf conf;
		conf.set_name(name);
		conf.set_type(engine + "_dropout");
		DropoutConf *dropout = conf.mutable_dropout_conf();
		dropout->set_dropout_ratio(dropout_ratio);

		return conf;
	}

	// Based on: https://github.com/chengyangfu/pytorch-vgg-cifar10/blob/master/vgg.py
	// // (also https://github.com/kuangliu/pytorch-cifar/blob/master/models/vgg.py)
	FeedForwardNet CreateNetVGG11()
	{
		// make layers
		FeedForwardNet net;
		Shape s{3, 32, 32};

		net.Add(GenConvConf("conv1", 64, 3, 1, 1), &s);
		net.Add(GenReLUConf("relu1"));
		net.Add(GenPoolingConf("pool1", true, 2, 2, 0));
		net.Add(GenConvConf("conv2", 128, 3, 1, 1));
		net.Add(GenReLUConf("relu2"));
		net.Add(GenPoolingConf("pool2", true, 2, 2, 0));
		net.Add(GenConvConf("conv3", 256, 3, 1, 1));
		net.Add(GenReLUConf("relu3"));
		net.Add(GenConvConf("conv4", 256, 3, 1, 1));
		net.Add(GenReLUConf("relu4"));
		net.Add(GenPoolingConf("pool3", true, 2, 2, 0));
		net.Add(GenConvConf("conv5", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu5"));
		net.Add(GenConvConf("conv6", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu6"));
		net.Add(GenPoolingConf("pool4", true, 2, 2, 0));
		net.Add(GenConvConf("conv7", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu7"));
		net.Add(GenConvConf("conv8", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu8"));
		net.Add(GenPoolingConf("pool5", true, 2, 2, 0));
		net.Add(GenFlattenConf("flat"));
		net.Add(GenDropoutConf("drop1", 0.5));
		net.Add(GenDenseConf("lin1", 512));
		net.Add(GenReLUConf("relu9"));
		net.Add(GenDropoutConf("drop2", 0.5));
		net.Add(GenDenseConf("lin2", 512));
		net.Add(GenReLUConf("relu10"));
		net.Add(GenDenseConf("lin3", 10));
		return net;
	}

	FeedForwardNet CreateNetVGG16()
	{
		// make layers
		FeedForwardNet net;
		Shape s{3, 32, 32};

		net.Add(GenConvConf("conv1", 64, 3, 1, 1), &s);
		net.Add(GenReLUConf("relu1"));
		net.Add(GenConvConf("conv2", 64, 3, 1, 1));
		net.Add(GenReLUConf("relu2"));
		net.Add(GenPoolingConf("pool1", true, 2, 2, 0));
		net.Add(GenConvConf("conv3", 128, 3, 1, 1));
		net.Add(GenReLUConf("relu3"));
		net.Add(GenConvConf("conv4", 128, 3, 1, 1));
		net.Add(GenReLUConf("relu4"));
		net.Add(GenPoolingConf("pool2", true, 2, 2, 0));
		net.Add(GenConvConf("conv5", 256, 3, 1, 1));
		net.Add(GenReLUConf("relu5"));
		net.Add(GenConvConf("conv6", 256, 3, 1, 1));
		net.Add(GenReLUConf("relu6"));
		net.Add(GenConvConf("conv7", 256, 3, 1, 1));
		net.Add(GenReLUConf("relu7"));
		net.Add(GenPoolingConf("pool3", true, 2, 2, 0));
		net.Add(GenConvConf("conv8", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu8"));
		net.Add(GenConvConf("conv9", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu9"));
		net.Add(GenConvConf("conv10", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu10"));
		net.Add(GenPoolingConf("pool4", true, 2, 2, 0));
		net.Add(GenConvConf("conv11", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu11"));
		net.Add(GenConvConf("conv12", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu12"));
		net.Add(GenConvConf("conv13", 512, 3, 1, 1));
		net.Add(GenReLUConf("relu13"));		
		net.Add(GenPoolingConf("pool5", true, 2, 2, 0));
		net.Add(GenFlattenConf("flat"));
		net.Add(GenDropoutConf("drop1", 0.5));
		net.Add(GenDenseConf("lin1", 512));
		net.Add(GenReLUConf("relu9"));
		net.Add(GenDropoutConf("drop2", 0.5));
		net.Add(GenDenseConf("lin2", 512));
		net.Add(GenReLUConf("relu10"));
		net.Add(GenDenseConf("lin3", 10));
		return net;
	}

	void Train(int num_epoch, string data_dir)
	{
		LOG(INFO) << "*** Train";

		Cifar10 data(data_dir);
		Tensor train_x, train_y, test_x, test_y;
		{
			auto train = data.ReadTrainData();
			size_t nsamples = train.first.shape(0);
			auto mtrain =
				Reshape(train.first, Shape{nsamples, train.first.Size() / nsamples});
			const Tensor &mean = Average(mtrain, 0);
			SubRow(mean, &mtrain);
			train_x = Reshape(mtrain, train.first.shape());
			train_y = train.second;
			auto test = data.ReadTestData();
			nsamples = test.first.shape(0);
			auto mtest =
				Reshape(test.first, Shape{nsamples, test.first.Size() / nsamples});
			SubRow(mean, &mtest);
			test_x = Reshape(mtest, test.first.shape());
			test_y = test.second;
		}
		CHECK_EQ(train_x.shape(0), train_y.shape(0));
		CHECK_EQ(test_x.shape(0), test_y.shape(0));
		LOG(INFO) << "Training samples = " << train_y.shape(0)
				  << ", Test samples = " << test_y.shape(0);
		auto net = CreateNetVGG16();
		SGD sgd;
		OptimizerConf opt_conf;
		opt_conf.set_momentum(0.9);
		auto reg = opt_conf.mutable_regularizer();
		reg->set_coefficient(0.0005);
		sgd.Setup(opt_conf);
		sgd.SetLearningRateGenerator([](int step) {
			return 0.01f / static_cast<float>(1u << (step/30));
		});

		SoftmaxCrossEntropy loss;
		Accuracy acc;
		net.Compile(true, &sgd, &loss, &acc);
#ifdef USE_CUDNN
		auto dev = std::make_shared<CudaGPU>();
		net.ToDevice(dev);
		train_x.ToDevice(dev);
		train_y.ToDevice(dev);
		test_x.ToDevice(dev);
		test_y.ToDevice(dev);
#endif // USE_CUDNN
		net.Train(100, num_epoch, train_x, train_y, test_x, test_y);
		test_y = net.Forward(kEval, test_x);

#ifdef MY_FILE_READER
		FileReader snap;
		snap.OpenForWrite("myfilesnap.bin");
#else
		Snapshot snap("mysnap", Snapshot::kWrite, 100);
#endif
		vector<string> names = net.GetParamNames();
		vector<Tensor> params = net.GetParamValues();
		vector<string>::iterator namesIter;
		vector<Tensor>::iterator paramsIter;
		int idx = 0;
		for (paramsIter = params.begin(), namesIter = names.begin();
			 paramsIter != params.end() && namesIter != names.end();
			 paramsIter++, namesIter++)
		{
			LOG(INFO) << "Write param: " << *namesIter;
			snap.Write(*namesIter, *paramsIter);
			idx++;
		}
#ifdef MY_FILE_READER
		snap.Close();
#endif
	}

	void Eval(string data_dir)
	{
		LOG(INFO) << "*** Evaluate";

		Cifar10 data(data_dir);
		Tensor test_x, test_y;
		{
			auto test = data.ReadTestData();
			size_t nsamples = test.first.shape(0);
			auto mtest =
				Reshape(test.first, Shape{nsamples, test.first.Size() / nsamples});
			const Tensor &mean = Average(mtest, 0);
			SubRow(mean, &mtest);
			test_x = Reshape(mtest, test.first.shape());
			test_y = test.second;
		}
		CHECK_EQ(test_x.shape(0), test_y.shape(0));
		LOG(INFO) << "Test samples = " << test_y.shape(0);
		auto net = CreateNetVGG16();
		SGD sgd;
		OptimizerConf opt_conf;
		opt_conf.set_momentum(0.9);
		auto reg = opt_conf.mutable_regularizer();
		reg->set_coefficient(0.004);
		sgd.Setup(opt_conf);
		sgd.SetLearningRateGenerator([](int step)
									 {
										 if (step <= 120)
											 return 0.001;
										 else if (step <= 130)
											 return 0.0001;
										 else
											 return 0.00001;
									 });

		SoftmaxCrossEntropy loss;
		Accuracy acc;
		net.Compile(true, &sgd, &loss, &acc);

#ifdef MY_FILE_READER
		FileReader snap;
		snap.OpenForRead("myfilesnap.bin");
#else
		Snapshot snap("mysnap", Snapshot::kRead, 100);
#endif
		vector<std::pair<std::string, Tensor>> params = snap.Read();
		net.SetParamValues(params);

#ifdef USE_CUDNN
		auto dev = std::make_shared<CudaGPU>();
		net.ToDevice(dev);
		test_x.ToDevice(dev);
		test_y.ToDevice(dev);
#endif // USE_CUDNN

		float val = 0.f;
		LOG(INFO) << "Start evaluation";
		std::pair<Tensor, Tensor> ret = net.EvaluateOnBatchAccuracy(test_x, test_y, &val);
		LOG(INFO) << "Accuracy: " << val;
#ifdef MY_FILE_READER
		snap.Close();
#endif
	}
}

int main(int argc, char **argv)
{
	singa::InitChannel(nullptr);
	int pos = singa::ArgPos(argc, argv, "-epoch");
	int nEpoch = 1;
	if (pos != -1)
		nEpoch = atoi(argv[pos + 1]);
	pos = singa::ArgPos(argc, argv, "-data");
	string data = "cifar-10-batches-bin";
	if (pos != -1)
		data = argv[pos + 1];

	LOG(INFO) << "Start training";
	singa::Train(nEpoch, data);
	LOG(INFO) << "End training";

	LOG(INFO) << "Start evaluation";
	singa::Eval(data);
	LOG(INFO) << "End evaluation";
}
