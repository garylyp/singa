#include <string>

#include "singa/model/feed_forward_net.h"
#include "singa/model/optimizer.h"
#include "singa/model/metric.h"
#include "singa/utils/channel.h"
#include "singa/utils/string.h"
#include "singa/io/snapshot.h"

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
static const size_t kBatchSize = 1;
static const size_t kImageDim = 32;

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
	ReLUConf *relu = conf.mutable_relu_conf();	
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
	bspec->set_lr_mult(2);
	bspec->set_decay_mult(0);

	return conf;
}

LayerConf GenLRNConf(string name) {
	LayerConf conf;
	conf.set_name(name);
	conf.set_type(engine + "_lrn");
	LRNConf *lrn = conf.mutable_lrn_conf();
	lrn->set_local_size(3);
	lrn->set_alpha(5e-05);
	lrn->set_beta(0.75);
	return conf;
}

LayerConf GenFlattenConf(string name) {
	LayerConf conf;
	conf.set_name(name);
	conf.set_type("singa_flatten");
	FlattenConf *flat = conf.mutable_flatten_conf();
	return conf;
}

FeedForwardNet CreateNet() {
	FeedForwardNet net;
	Shape s{3, 32, 32};

	net.Add(GenConvConf("conv1", 32, 5, 1, 2, 0.0001), &s);
	net.Add(GenPoolingConf("pool1", true, 3, 2, 1));
	net.Add(GenReLUConf("relu1"));
	net.Add(GenLRNConf("lrn1"));
	net.Add(GenConvConf("conv2", 32, 5, 1, 2, 0.01));
	net.Add(GenReLUConf("relu2"));
	net.Add(GenPoolingConf("pool2", false, 3, 2, 1));
	net.Add(GenLRNConf("lrn2"));
	net.Add(GenConvConf("conv3", 64, 5, 1, 2, 0.01));
	net.Add(GenReLUConf("relu3"));
	net.Add(GenPoolingConf("pool3", false, 3, 2, 1));
	net.Add(GenFlattenConf("flat"));
	net.Add(GenDenseConf("ip", 10, 0.01, 250));
	return net;
}
/* TODO remove
const std::pair<Tensor, Tensor> ReadImageFile(string file) {	
	Tensor image(Shape{1, 3, kImageDim, kImageDim});
	Tensor label(Shape{1}, kInt);
	std::ifstream data_file(file.c_str(), std::ios::in | std::ios::binary);
	if (!data_file.is_open())
		LOG(ERROR) << "Unable to open file " << file;
	LOG(INFO) << "Read image from file";
	char buff[kImageVol + 1];
	float data[kImageVol];
	data_file.read(buff, kImageVol + 1);
	int label_val = (int)buff[0];
	label.CopyDataFromHostPtr(&label_val, 1, 0);
	for (int i = 0; i < kImageVol; i++)
		data[i] = (float)((int)buff[i+1]);
	image.CopyDataFromHostPtr(data, kImageVol, 0);
	data_file.close();
	return std::make_pair(image, label);
}
*/

std::pair<Tensor, Tensor> LoadData() {
	LOG(INFO) << "Load image from memory";
	Tensor images(Shape{1, 3, kImageSize, kImageSize});
	Tensor labels(Shape{1}, kInt);
	char* ptr = (char*)(&_binary_images_bin_start);
	float data[kImageVol];
	for (int i = 0; i < kImageVol; i++)
		data[i] = (float)((int)ptr[i+1]);
	images.CopyDataFromHostPtr(data, kImageVol, 0);	
	int lbl = (int)(*ptr);
	labels.CopyDataFromHostPtr(&lbl, 1, 0);
	return std::make_pair(images, labels);
}

vector<std::pair<std::string, Tensor>> LoadParams() {
	std::unordered_set<std::string> param_names_;
	std::unordered_map<std::string, Tensor> param_map_;
	std::string key;
	char* ptr;
	size_t size; 

	int param_size = &_binary_myfilesnap_bin_end - &_binary_myfilesnap_bin_start;
	LOG(INFO) << "Size of parameters: " << param_size;

	MemReader mem_reader((char*)&_binary_myfilesnap_bin_start, param_size);

	return mem_reader.Read();

	while (mem_reader.Read(&key, &ptr, &size)) {
		CHECK(param_names_.count(key) == 0);
		LOG(INFO) << "Read param: " << key;
		param_names_.insert(key);		
		param_map_[key].FromBytes((uint8_t *)ptr, size);
	}
	std::vector<std::pair<std::string, Tensor>> ret;
	for (auto it = param_map_.begin(); it != param_map_.end(); ++it)
		ret.push_back(*it);
	return ret;
}

int Eval() {
	Tensor test_x, test_y;
	
	auto test = LoadData();
	// auto test = ReadImageFile("images.bin");
	test_x = test.first;
	test_y = test.second;		
	CHECK_EQ(test_x.shape(0), test_y.shape(0));
	LOG(INFO) << "Test samples = " << test_y.shape(0);

	SGD sgd;
	OptimizerConf opt_conf;
	opt_conf.set_momentum(0.9);
	auto reg = opt_conf.mutable_regularizer();
	reg->set_coefficient(0.004);
	sgd.Setup(opt_conf);
	sgd.SetLearningRateGenerator([](int step) {
		if (step <= 120)
			return 0.001;
		else if (step <= 130)
			return 0.0001;
		else
			return 0.00001;
	});

	LOG(INFO) << "Create net.";
	SoftmaxCrossEntropy loss;
	Accuracy acc;
	auto net = CreateNet();
	LOG(INFO) << "Compile net.";
	net.Compile(true, &sgd, &loss, &acc);

	LOG(INFO) << "Done create net.";

	vector<std::pair<std::string, Tensor>> params = LoadParams();
	net.SetParamValues(params);
	Tensor tout =  net.EvaluateOnBatchOutput(test_x, test_y);

	float vals[10];
	tout.GetValue(vals, 10);
	float max = vals[0];
	int max_pos = 0;
	LOG(INFO) << 0 << " " << vals[0];
	for (int i = 1; i < 10; i++) {
		if (max < vals[i]) {
			max = vals[i];
			max_pos = i;
		}
		LOG(INFO) << i << " " << vals[i];
	}
	return max_pos;
}

} // end singa

int main(int argc, char **argv) {
	singa::InitChannel(nullptr);

	int idx = singa::Eval();
	LOG(INFO) << "Label: " << idx;

	return 0;
}
