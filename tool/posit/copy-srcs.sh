#!/bin/bash

set -x

DST="build_bare_posit"

cd ../../
rm -rf $DST
mkdir $DST
cd $DST

cp ../tool/posit/Makefile .

cp ../examples/cifar10/cnn-one-image-bare.cc .
# cp ../examples/cifar10/cnn-bare.cc .
cp ../examples/cifar10/mem_reader.h .
cp ../examples/cifar10/mem_reader.cc .
cp ../examples/cifar10/objects.h .

cp ../examples/cifar10/images.o .
cp ../examples/cifar10/labels.o .
cp ../examples/cifar10/myfilesnap.o .

cp ../src/core/tensor/* .
cp ../src/core/device/device.cc .
cp ../src/core/device/platform.cc .
cp ../src/core/device/cpp_cpu.cc .
cp ../src/io/binfile* .
cp ../src/io/textfile_writer.cc .
cp ../src/model/layer/* .
cp ../src/model/optimizer/* .
cp ../src/model/loss/* .
cp ../src/model/metric/* .
cp ../src/model/updater/* .
cp ../src/utils/* .
cp ../src/model/feed_forward_net.cc .

cp -r ../include/singa .
rsync -a --ignore-existing ../build/include/singa/ singa/

rm -f cudnn*
rm -f opencl*
rm -f *.cu
rm -f *.cu
rm -f *.cl
rm -f tensor_math_opencl.h tensor_math_cuda.h
rm -f prelu.cc prelu.h rnn.cc rnn.h slice.cc slice.h dropout.cc dropout.h concat.cc concat.h \
batchnorm.cc batchnorm.h local_all_reduce.* local_updater.*
rm -rf singa/proto
