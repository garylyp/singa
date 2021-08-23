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

#include "singa/model/feed_forward_net.h"
#include "singa/model/optimizer.h"
#include "singa/model/metric.h"
#include "singa/utils/channel.h"
#include "singa/utils/string.h"
#include "singa/io/snapshot.h"
#include "singa/core/common.h"

#include "mem_reader.h"
#include "objects.h"

namespace singa {
// currently supports 'cudnn' and 'singacpp'
#ifdef USE_CUDNN
const std::string engine = "cudnn";
#else
const std::string engine = "singacpp";
#endif  // USE_CUDNN

static const size_t kImageSize = 32;
static const size_t kImageVol = 3072;
// batch size is 10000 but we may need to lower it to run on Freedom
static const size_t kBatchSize = 1000;

static const float my_alpha = 5e-05;
static const float my_beta = 0.75;
static const float my_momentum = 0.9;
static const float my_reg_coef = 0.004;

static const float my_10_minus_2 = 0.01;
static const float my_10_minus_3 = 0.001;
static const float my_10_minus_4 = 0.0001;
static const float my_10_minus_5 = 0.00001;


LayerConf GenConvConf(string name, int nb_filter, int kernel, int stride,
		int pad, float std) {
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
	wfill->set_std(std);

	ParamSpec *bspec = conf.add_param();
	bspec->set_name(name + "_bias");
	bspec->set_lr_mult(2);
	//  bspec->set_decay_mult(0);
	return conf;
}

LayerConf GenPoolingConf(string name, bool max_pool, int kernel, int stride,
		int pad) {
	LayerConf conf;
	conf.set_name(name);
	conf.set_type(engine + "_pooling");
	PoolingConf *pool = conf.mutable_pooling_conf();
	pool->set_kernel_size(kernel);
	pool->set_stride(stride);
	pool->set_pad(pad);
	if (!max_pool) pool->set_pool(PoolingConf_PoolMethod_AVE);
	return conf;
}

LayerConf GenReLUConf(string name) {
	LayerConf conf;
	conf.set_name(name);
	conf.set_type(engine + "_relu");
	conf.mutable_relu_conf();
	return conf;
}

LayerConf GenDenseConf(string name, int num_output, float std, float wd) {
	LayerConf conf;
	conf.set_name(name);
	conf.set_type("singa_dense");
	DenseConf *dense = conf.mutable_dense_conf();
	dense->set_num_output(num_output);

	ParamSpec *wspec = conf.add_param();
	wspec->set_name(name + "_weight");
	wspec->set_decay_mult(wd);
	auto wfill = wspec->mutable_filler();
	wfill->set_type("Gaussian");
	wfill->set_std(std);

	ParamSpec *bspec = conf.add_param();
	bspec->set_name(name + "_bias");
	bspec->set_lr_mult(const_float_two);
	bspec->set_decay_mult(const_float_zero);

	return conf;
}

LayerConf GenLRNConf(string name) {
	LayerConf conf;
	conf.set_name(name);
	conf.set_type(engine + "_lrn");
	LRNConf *lrn = conf.mutable_lrn_conf();
	lrn->set_local_size(3);
	lrn->set_alpha(my_alpha);
	lrn->set_beta(my_beta);
	return conf;
}

LayerConf GenFlattenConf(string name) {
	LayerConf conf;
	conf.set_name(name);
	conf.set_type("singa_flatten");
	conf.mutable_flatten_conf();
	return conf;
}

FeedForwardNet CreateNet() {
	FeedForwardNet net;
	Shape s{3, 32, 32};

	net.Add(GenConvConf("conv1", 32, 5, 1, 2, my_10_minus_4), &s);
	net.Add(GenPoolingConf("pool1", true, 3, 2, 1));
	net.Add(GenReLUConf("relu1"));
	net.Add(GenLRNConf("lrn1"));
	net.Add(GenConvConf("conv2", 32, 5, 1, 2, my_10_minus_2));
	net.Add(GenReLUConf("relu2"));
	net.Add(GenPoolingConf("pool2", false, 3, 2, 1));
	net.Add(GenLRNConf("lrn2"));
	net.Add(GenConvConf("conv3", 64, 5, 1, 2, my_10_minus_2));
	net.Add(GenReLUConf("relu3"));
	net.Add(GenPoolingConf("pool3", false, 3, 2, 1));
	net.Add(GenFlattenConf("flat"));
	net.Add(GenDenseConf("ip", 10, my_10_minus_2, 250));
	return net;
}

std::pair<Tensor, Tensor> LoadData() {
	Tensor images(Shape{kBatchSize, 3, kImageSize, kImageSize});
	Tensor labels(Shape{kBatchSize}, kInt);
	for (size_t itemid = 0; itemid < kBatchSize; ++itemid) {
		images.CopyDataFromHostPtr((float*)&_binary_images_bin_start + (itemid * kImageVol), kImageVol, itemid * kImageVol);
	}
	labels.CopyDataFromHostPtr((int*)(&_binary_labels_bin_start), kBatchSize);
	return std::make_pair(images, labels);
}

vector<std::pair<std::string, Tensor>> LoadParams() {
#ifndef LITE_POSIT
	std::unordered_set<std::string> param_names_;
	std::unordered_map<std::string, Tensor> param_map_;
	singa::TensorProto tp;
	std::string key, val;

	int param_size = &_binary_mysnap_bin_end - &_binary_mysnap_bin_start;
	LOG(INFO) << "Size of parameters: " << param_size;

	MemReader mem_reader((char*)&_binary_mysnap_bin_start, param_size);

	return mem_reader.Read();

	while (mem_reader.Read(&key, &val)) {
		CHECK(param_names_.count(key) == 0);
		LOG(INFO) << "Read param: " << key;
		param_names_.insert(key);
		CHECK(tp.ParseFromString(val));
		param_map_[key].FromProto(tp);
	}
	std::vector<std::pair<std::string, Tensor>> ret;
	for (auto it = param_map_.begin(); it != param_map_.end(); ++it)
		ret.push_back(*it);
	return ret;
#else
	int param_size = &_binary_myfilesnap_bin_end - &_binary_myfilesnap_bin_start;
	LOG(INFO) << "Size of parameters: " << param_size;
	MemReader mem_reader((char*)&_binary_myfilesnap_bin_start, param_size);
	return mem_reader.Read();
#endif
}

void Eval() {
	Tensor test_x, test_y;
	{
		auto test = LoadData();
		size_t nsamples = test.first.shape(0);
		auto mtest =
				Reshape(test.first, Shape{nsamples, test.first.Size() / nsamples});
		const Tensor& mean = Average(mtest, 0);
		SubRow(mean, &mtest);
		test_x = Reshape(mtest, test.first.shape());
		test_y = test.second;
	}
	CHECK_EQ(test_x.shape(0), test_y.shape(0));
	LOG(INFO) << "Test samples = " << test_y.shape(0);
	auto net = CreateNet();
	SGD sgd;
	OptimizerConf opt_conf;
	opt_conf.set_momentum(my_momentum);
	auto reg = opt_conf.mutable_regularizer();
	reg->set_coefficient(my_reg_coef);
	sgd.Setup(opt_conf);
	sgd.SetLearningRateGenerator([](int step) {
		if (step <= 120)
			return my_10_minus_3;
		else if (step <= 130)
			return my_10_minus_4;
		else
			return my_10_minus_5;
	});

	SoftmaxCrossEntropy loss;
	Accuracy acc;
	net.Compile(true, &sgd, &loss, &acc);

	vector<std::pair<std::string, Tensor>> params = LoadParams();
	net.SetParamValues(params);
	float val = const_float_zero;
	std::pair<Tensor, Tensor> ret = net.EvaluateOnBatchAccuracy(test_x, test_y, &val);
	LOG(INFO) << "Accuracy: " << val;
}

}

int main(int argc, char **argv) {
	singa::InitChannel(nullptr);

	LOG(INFO) << "Start evaluation";
	singa::Eval();
	LOG(INFO) << "End evaluation";
}
