/**
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SINGA_MODEL_OPTIMIZER_CONF__H_
#define SINGA_MODEL_OPTIMIZER_CONF__H_

#include <string>
#include <vector>
#include "singa/core/common.h"
#include "singa/utils/logging.h"

namespace singa {

enum Phase : int {
	kTrain = 4,
			kEval = 8
};

enum PoolingConf_PoolMethod : int {
	MAX = 0,
			PoolingConf_PoolMethod_MAX = 0,
			AVE = 1,
			PoolingConf_PoolMethod_AVE = 1,
			STOCHASTIC = 2,
			PoolingConf_PoolMethod_STOCHASTIC =2
};

class ConstraintConf {
public:
	~ConstraintConf();

	bool ParseFromString(const std::string& data) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}

	inline std::string type() const {
		return type_;
	}

	inline void set_type(std::string value) {
		type_ = value;
	}

	inline float threshold() const {
		return threshold_;
	}

	inline void set_threshold(float value) {
		threshold_ = value;
	}


private:
	// case insensitive to limit the parameter value/gradient scale
	std::string type_ = "l2";
	// e.g., the threshold for limiting the parameter scale.
	float threshold_ = 0.0;

	void CopyFrom(const ConstraintConf& from);

};

class RegularizerConf {
public:
	~RegularizerConf();

	bool ParseFromString(const std::string& data) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}

	inline std::string type() const {
		return type_;
	}

	inline void set_tyep(std::string value) {
		type_ = value;
	}

	inline float coefficient() const {
		return coefficient_;
	}

	inline void set_coefficient(float value) {
		coefficient_ = value;
	}

private:

	// case insensitive to regularize the parameters, e.g., L2.
	std::string type_ = "l2";
	// e.g., the weight decay for L2 regularizer
	float coefficient_ = 0.0;

	void CopyFrom(const RegularizerConf& from);
};

class OptimizerConf {
public:
	~OptimizerConf();

	bool ParseFromString(const std::string& data) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}

	inline RegularizerConf* mutable_regularizer() {
		if (regularizer_ == nullptr) {
			auto* p = new RegularizerConf();
			regularizer_ = p;
		}
		return regularizer_;
	}

	inline RegularizerConf regularizer() const {
		CHECK_NOTNULL(regularizer_);
		return *regularizer_;
	}

	inline bool has_regularizer() const {
		return (regularizer_ != nullptr);
	}

	inline ConstraintConf constraint() const {
		CHECK_NOTNULL(constraint_);
		return *constraint_;
	}

	inline bool has_constraint() const {
		return (constraint_ != nullptr);
	}

	inline bool has_momentum() const {
		return true;
	}

	inline float momentum() const {
		return momentum_;
	}

	inline void set_momentum(float value) {
		momentum_ = value;
	}

	inline float rho() const {
		return rho_;
	}

	inline float delta() const {
		return delta_;
	}

private:
	// case insensitive
	std::string type_ = "sgd";

	// used by RMSprop and Adadelta
	float rho_ = 0.95;

	// used by Adam and AdamMax
	float beta_1_ = 0.9;
	float beta_2_ = 0.999;

	// used by vanilla sgd and nesterov
	float momentum_ = 0.9;

	// delta is used to avoid dividing zero
	float delta_ = 1e-8;

	// global regularizer lower priority than ParamSpec regularizer
	RegularizerConf *regularizer_ = nullptr;
	// global constraint lower priority than ParamSpec constraint
	ConstraintConf *constraint_ = nullptr;

	void CopyFrom(const OptimizerConf& from);
};

class LossConf {
public:
	bool ParseFromString(const std::string& data) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}
};

class FillerConf {
public:
	bool ParseFromString(const std::string& data) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}

	inline void set_type(std::string value) {
		type_ = value;
	}

	inline std::string type() const {
		return type_;
	}

	inline void set_std(float value) {
		std_ = value;
	}

	inline float std() const {
		return std_;
	}

	inline float value() const {
		return value_;
	}

	inline float min() const {
		return min_;
	}

	inline float max() const {
		return max_;
	}

	inline float mean() const {
		return mean_;
	}

private:
	// The filler type, case insensitive
	std::string type_ = "constant";
	float value_ = const_float_zero; // the value in constant filler
	float min_ = const_float_zero; // the min value in uniform filler
	float max_ = const_float_one; // the max value in uniform filler
	float mean_ = const_float_zero; // the mean value in Gaussian filler
	float std_ = const_float_one; // the std value in Gaussian filler
	// The expected number of non-zero output weights for a given input in
	// Gaussian filler -- the default -1 means don't perform sparsification.
	/* optional int32 sparse = 7 [default = -1]; */
	// Normalize the filler variance by fan_in, fan_out, or their average.
	// Applies to 'xavier' and 'msra' fillers.
	enum VarianceNorm {
		FAN_IN = 0,
		FAN_OUT = 1,
		AVERAGE = 2
	};
	VarianceNorm variance_norm_ = FAN_IN;
};

class ParamSpec {
public:
	ParamSpec();
	~ParamSpec();

	bool ParseFromString(const std::string& data) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}

	inline std::string name() const {
		return name_;
	}

	inline void set_name(std::string value) {
		name_ = value;
	}

	inline bool has_constraint() const {
		return (constraint_ != nullptr);
	}

	inline ConstraintConf constraint() const {
		CHECK_NOTNULL(constraint_);
		return *constraint_;
	}

	inline bool has_regularizer() const {
		return (regularizer_ != nullptr);
	}

	inline RegularizerConf regularizer() const {
		CHECK_NOTNULL(regularizer_);
		return *regularizer_;
	}

	inline FillerConf filler() const {
		CHECK_NOTNULL(filler_);
		return *filler_;
	}

	inline bool has_filler() const {
		return (filler_ != nullptr);
	}

	inline bool has_decay_mult() const {
		return true;
	}

	inline float decay_mult() const {
		return decay_mult_;
	}

	inline bool has_lr_mult() const {
		return true;
	}

	inline float lr_mult() const {
		return lr_mult_;
	}

	inline void set_lr_mult(float value) {
		lr_mult_ = value;
	}

	inline void set_decay_mult(float value) {
		decay_mult_ = value;
	}

	inline RegularizerConf* mutable_regularizer() {
		if (regularizer_ == nullptr)
			regularizer_ = new RegularizerConf();
		return regularizer_;
	}

	inline FillerConf* mutable_filler() {
		if (filler_ == nullptr)
			filler_ = new FillerConf();
		return filler_;
	}

private:
	// The names of the parameter blobs -- useful for sharing parameters among
	// layers, but never required otherwise.  To share a parameter between two
	// layers, give it a (non-empty) name.
	std::string name_;

	// The multiplier on the global learning rate for this parameter.
	float lr_mult_ = 1.0;

	// The multiplier on the global weight decay for this parameter.
	float decay_mult_ = 1.0;

	// SINGA uses this filed internally. Users just configure the fillers in
	// Layer specific conf message as caffe (style).
	FillerConf* filler_ = nullptr;
	ConstraintConf* constraint_ = nullptr;
	RegularizerConf* regularizer_ = nullptr;

	void CopyFrom(const ParamSpec& from);
};

class BlobProto {

};

class DenseConf {
public:
	inline uint32_t num_output() const {
		return num_output_;
	}

	inline bool transpose() const {
		return transpose_;
	}

	inline bool bias_term() const {
		return bias_term_;
	}

	inline void set_num_output(uint32_t value) {
		num_output_ = value;
	}

private:
	uint32_t num_output_ = 0; // The number of outputs for the layer
	bool bias_term_ = true; // whether to have bias terms
	FillerConf* weight_filler_ = nullptr; // The filler for the weight
	FillerConf* bias_filler_ = nullptr; // The filler for the bias

	// The first axis to be lumped into a single inner product computation;
	// all preceding axes are retained in the output.
	// May be negative to index from the end (e.g., -1 for the last axis).
	int32_t axis_ = 1;

	bool transpose_ = false; // whether transpose or not
};

class MetricConf {
public:

	bool ParseFromString(const std::string& data) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}

	inline uint32_t top_k() const {
		return top_k_;
	}

private:
	// When computing accuracy, count as correct by comparing the true label to
	// the top k scoring classes.  By default, only compare to the top scoring
	// class (i.e. argmax).
	uint32_t top_k_ = 1;

	// The "label" axis of the prediction blob, whose argmax corresponds to the
	// predicted label -- may be negative to index from the end (e.g., -1 for the
	// last axis).  For example, if axis == 1 and the predictions are
	// (N x C x H x W), the label blob is expected to contain N*H*W ground truth
	// labels with integer values in {0, 1, ..., C-1}.
	int32_t axis = 1;

	// If specified, ignore instances with the given label.
	int32_t ignore_label_ = 3;
};

class BatchNormConf {

};

class SplitConf {
public :
	inline int32_t output_size() const {
		return output_size_;
	}

private:
	int32_t output_size_ = 2;
};

class ReLUConf {
public:

	inline float negative_slope() const {
		return negative_slope_;
	}

private:
	// Allow non-zero slope for negative inputs to speed up optimization
	// Described in:
	// Maas, A. L., Hannun, A. Y., & Ng, A. Y. (2013). Rectifier nonlinearities
	// improve neural network acoustic models. In ICML Workshop on Deep Learning
	// for Audio, Speech, and Language Processing.
	float negative_slope_ = 0.0;

};

class ConvolutionConf {
public:

	inline size_t kernel_size_size() const {
		return kernel_size_.size();
	}

	inline size_t pad_size() const {
		return pad_.size();
	}

	inline size_t stride_size() const {
		return stride_.size();
	}

	inline uint32_t num_output() const {
		return num_output_;
	}

	inline void set_num_output(uint32_t value) {
		num_output_ = value;
	}

	inline uint32_t kernel_size(size_t i) const {
		CHECK_LT(i, kernel_size_.size());
		return kernel_size_.at(i);
	}

	inline void add_kernel_size(uint32_t value) {
		kernel_size_.push_back(value);
	}

	inline uint32_t pad(size_t i) const {
		CHECK_LT(i, pad_.size());
		return pad_.at(i);
	}

	inline void add_pad(uint32_t value) {
		pad_.push_back(value);
	}

	inline uint32_t stride(size_t i) const {
		CHECK_LT(i, stride_.size());
		return stride_.at(i);
	}

	inline void add_stride(uint32_t value) {
		stride_.push_back(value);
	}

	inline uint32_t kernel_w() const {
		return kernel_w_;
	}

	inline uint32_t kernel_h() const {
		return kernel_h_;
	}

	inline uint32_t pad_w() const {
		return pad_w_;
	}

	inline uint32_t pad_h() const {
		return pad_h_;
	}

	inline uint32_t stride_w() const {
		return stride_w_;
	}

	inline uint32_t stride_h() const {
		return stride_h_;
	}

	inline bool has_stride_w() const {
		return (stride_w_ != 0);
	}

	inline bool has_stride_h() const {
		return (stride_h_ != 0);
	}

	inline bool bias_term() const {
		return bias_term_;
	}

	inline void set_bias_term(bool value) {
		bias_term_ = value;
	}

private:
	uint32_t num_output_ = 0; // The number of outputs for the layer
	bool bias_term_ = true; // whether to have bias terms

	// Pad, kernel size, and stride are all given as a single value for equal
	// dimensions in all spatial dimensions, or once per spatial dimension.
	std::vector<uint32_t> pad_; // The padding size; defaults to 0
	std::vector<uint32_t> kernel_size_; // The kernel size
	std::vector<uint32_t> stride_; // The stride; defaults to 1

	// For 2D convolution only, the *_h and *_w versions may also be used to
	// specify both spatial dimensions.
	uint32_t pad_h_ = 0; // The padding height (2D only)
	uint32_t pad_w_ = 0; // The padding width (2D only)
	uint32_t kernel_h_ = 0; // The kernel height (2D only)
	uint32_t kernel_w_ = 0; // The kernel width (2D only)
	uint32_t stride_h_ = 0; // The stride height (2D only)
	uint32_t stride_w_ = 0; // The stride width (2D only)

	// SINGA: not supported.
	// optional uint32 group = 5 [default = 1]; // The group size for group conv

	FillerConf* weight_filler_ = nullptr; // The filler for the weight
	FillerConf* bias_filler_ = nullptr; // The filler for the bias
	enum Engine {
		DEFAULT = 0,
		CAFFE = 1,
		CUDNN = 2
	};
	Engine engine = DEFAULT;

	// The axis to interpret as "channels" when performing convolution.
	// Preceding dimensions are treated as independent inputs;
	// succeeding dimensions are treated as "spatial".
	// With (N, C, H, W) inputs, and axis == 1 (the default), we perform
	// N independent 2D convolutions, sliding C-channel (or (C/g)-channels, for
	// groups g>1) filters across the spatial axes (H, W) of the input.
	// With (N, C, D, H, W) inputs, and axis == 1, we perform
	// N independent 3D convolutions, sliding (C/g)-channels
	// filters across the spatial axes (D, H, W) of the input.
	// SINGA: not supported;
	// optional int32 axis = 16 [default = 1];

	// Whether to force use of the general ND convolution, even if a specific
	// implementation for blobs of the appropriate number of spatial dimensions
	// is available. (Currently, there is only a 2D-specific convolution
	// implementation; for input blobs with num_axes != 2, this option is
	// ignored and the ND implementation will be used.)
	// SINGA: not supported;
	// optional bool force_nd_im2col = 17 [default = false];


	// SINGA: add by xiangrui
	// cudnn workspace size in MB
	int32_t workspace_byte_limit_ = 8192;

	// cudnn algorithm preference
	// options: "fastest", "limited_workspace", "no_workspace"
	std::string prefer_ = "fastest";
};

class PoolingConf {
public:

	inline bool has_kernel_size() const {
		return true;
	}

	inline uint32_t kernel_size() const {
		return kernel_size_;
	}

	inline void set_kernel_size(uint32_t value) {
		kernel_size_ = value;
	}

	inline uint32_t kernel_h() const {
		return kernel_h_;
	}

	inline uint32_t kernel_w() const {
		return kernel_w_;
	}

	inline bool has_pad() const {
		return (pad_ != 0);
	}

	inline uint32_t pad() const {
		return pad_;
	}

	inline void set_pad(uint32_t value) {
		pad_ = value;
	}

	inline uint32_t pad_w() const {
		return pad_w_;
	}

	inline uint32_t pad_h() const {
		return pad_h_;
	}

	inline bool has_stride() const {
		return (stride_ != 0);
	}

	inline uint32_t stride() const {
		return stride_;
	}

	inline void set_stride(uint32_t value) {
		stride_ = value;
	}

	inline uint32_t stride_w() const {
		return stride_w_;
	}

	inline uint32_t stride_h() const {
		return stride_h_;
	}

	inline PoolingConf_PoolMethod pool() const {
		return pool_;
	}

	inline void set_pool (PoolingConf_PoolMethod value) {
		pool_ = value;
	}

	inline bool ceil() const {
		return ceil_;
	}

private:
	PoolingConf_PoolMethod pool_ =  MAX; // The pooling method
	// Pad, kernel size, and stride are all given as a single value for equal
	// dimensions in height and width or as Y, X pairs.
	uint32_t pad_ = 0; // The padding size (equal in Y, X)
	uint32_t pad_h_ = 0; // The padding height
	uint32_t pad_w_ = 0; // The padding width
	uint32_t kernel_size_ = 0; // The kernel size (square)
	uint32_t kernel_h_ = 0; // The kernel height
	uint32_t kernel_w_ = 0; // The kernel width
	uint32_t stride_ = 1; // The stride (equal in Y, X)
	uint32_t stride_h_ = 0; // The stride height
	uint32_t stride_w_ = 0; // The stride width

	// If global_pooling then it will pool over the size of the bottom by doing
	// kernel_h = bottom->height and kernel_w = bottom->width
	bool global_pooling_ = false;
	// whether to propagate nan
	bool nan_prop_ = false;

	// Added by xiangrui, 18 Oct, 2016
	bool ceil_ = false;
};

class LRNConf {
public:

	inline void set_local_size(uint32_t value) {
		local_size_ = value;
	}

	inline void set_alpha(float value) {
		alpha_ = value;
	}

	inline void set_beta(float value) {
		beta_ = value;
	}

	inline void set_k(float value) {
		k_ = value;
	}

	inline uint32_t local_size() const {
		return local_size_;
	}

	inline float alpha() const {
		return alpha_;
	}

	inline float beta() const {
		return beta_;
	}

	inline float k() const {
		return k_;
	}

private:
	uint32_t local_size_ = 5;
	float alpha_ = const_float_one;
	float beta_ = 0.75;
	enum NormRegion {
		ACROSS_CHANNELS = 0,
		WITHIN_CHANNEL = 1
	};
	NormRegion norm_region_ = ACROSS_CHANNELS;
	float k_ = const_float_one;
};

class FlattenConf {
public:

	inline int32_t axis() const {
		return axis_;
	}

	inline int32_t end_axis() const {
		return end_axis_;
	}

private:
	// The first axis to flatten: all preceding axes are retained in the output.
	// May be negative to index from the end (e.g., -1 for the last axis).
	int32_t axis_ = 1;

	// The last axis to flatten: all following axes are retained in the output.
	// May be negative to index from the end (e.g., the default -1 for the last
	// axis).
	int32_t end_axis_ = -1;
};

class LayerConf {
public:

	bool ParseFromString(const std::string& data) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}

	bool SerializeToString(std::string* output) {
		LOG(FATAL) << "ParseFromString not implemented";
		return false;
	}

	inline std::string name() const {
		return name_;
	}

	inline void set_name(std::string value) {
		has_name_ = true;
		name_ = value;
	}

	inline bool has_name() const {
		return has_name_;
	}

	inline std::string type() const {
		return type_;
	}

	inline void set_type(std::string value) {
		type_ = value;
	}

	inline std::vector<ParamSpec*> param_ptr() const {
		return params_;
	}

	inline std::vector<ParamSpec> param() const {
		std::vector<ParamSpec> params;
		for (auto param_spec_ptr : params_)
			params.push_back(*param_spec_ptr);
		return params;
	}

	inline DenseConf dense_conf() const {
		CHECK_NOTNULL(dense_conf_);
		return *dense_conf_;
	}

	inline SplitConf split_conf() const {
		CHECK_NOTNULL(split_conf_);
		return *split_conf_;
	}

	inline ReLUConf relu_conf() const {
		CHECK_NOTNULL(relu_conf_);
		return *relu_conf_;
	}

	inline PoolingConf pooling_conf() const {
		CHECK_NOTNULL(pooling_conf_);
		return *pooling_conf_;
	}

	inline ConvolutionConf convolution_conf() const {
		CHECK_NOTNULL(convolution_conf_);
		return *convolution_conf_;
	}

	inline LRNConf lrn_conf() const {
		CHECK_NOTNULL(lrn_conf_);
		return *lrn_conf_;
	}

	inline FlattenConf flatten_conf() const {
		CHECK_NOTNULL(flatten_conf_);
		return *flatten_conf_;
	}

	inline ConvolutionConf* mutable_convolution_conf() {
		if (convolution_conf_ == nullptr)
			convolution_conf_ = new ConvolutionConf();
		return convolution_conf_;
	}

	inline LRNConf* mutable_lrn_conf() {
		if (lrn_conf_ == nullptr)
			lrn_conf_ = new LRNConf();
		return lrn_conf_;
	}

	inline PoolingConf* mutable_pooling_conf() {
		if (pooling_conf_ == nullptr)
			pooling_conf_ = new PoolingConf();
		return pooling_conf_;
	}

	inline DenseConf* mutable_dense_conf() {
		if (dense_conf_ == nullptr)
			dense_conf_ = new DenseConf();
		return dense_conf_;
	}

	inline ReLUConf* mutable_relu_conf() {
		if (relu_conf_ == nullptr)
			relu_conf_ = new ReLUConf();
		return relu_conf_;
	}

	inline FlattenConf* mutable_flatten_conf() {
		if (flatten_conf_ == nullptr)
			flatten_conf_ = new FlattenConf();
		return flatten_conf_;
	}

	inline ParamSpec* add_param() {
		ParamSpec* new_param_spec = new ParamSpec();
		params_.push_back(new_param_spec);
		return new_param_spec;
	}

private:
	bool has_name_ = false;
	std::string name_;
	std::string type_;
	std::vector<ParamSpec*> params_;

	// The blobs containing the numeric parameters of the layer.
	std::vector<BlobProto*> blobs_;

	DenseConf* dense_conf_= nullptr;
	MetricConf* metric_conf_ = nullptr;
	BatchNormConf* batchnorm_conf_ = nullptr;
	SplitConf* split_conf_= nullptr;
	ReLUConf* relu_conf_ = nullptr;
	PoolingConf* pooling_conf_ = nullptr;
	ConvolutionConf* convolution_conf_ = nullptr;
	LRNConf* lrn_conf_ = nullptr;
	FlattenConf* flatten_conf_ = nullptr;
};

} // namespace singa

#endif // SINGA_MODEL_OPTIMIZER_CONF__H_
