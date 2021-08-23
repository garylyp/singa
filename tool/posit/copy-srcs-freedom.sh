#!/bin/bash

SRC="/media/dumi/home/dumi/git/riscv/incubator-singa"
DST="/media/dumi/home/dumi/git/riscv/freedom-e-sdk-dloghin/software/singa-cifar10"

LD_EXE="/media/dumi/home/dumi/git/riscv/riscv-tools/bin/riscv64-unknown-elf-ld"

mkdir -p $DST
cd $DST
rm -r *

# cp $SRC/tool/posit/Makefile-fp32 Makefile
cp $SRC/tool/posit/Makefile Makefile

cp $SRC/examples/cifar10/cnn-bare.cc .
cp $SRC/examples/cifar10/mem_reader.h .
cp $SRC/examples/cifar10/mem_reader.cc .
cp $SRC/examples/cifar10/objects.h .

cp $SRC/examples/cifar10/images.bin .
cp $SRC/examples/cifar10/labels.bin .
cp $SRC/examples/cifar10/myfilesnap.bin .

$LD_EXE -r -b binary images.bin -o images.o
$LD_EXE -r -b binary labels.bin -o labels.o
$LD_EXE -r -b binary myfilesnap.bin -o myfilesnap.o

cp $SRC/src/core/tensor/* .
cp $SRC/src/core/device/device.cc .
cp $SRC/src/core/device/platform.cc .
cp $SRC/src/core/device/cpp_cpu.cc .
cp $SRC/src/io/binfile* .
cp $SRC/src/io/textfile_writer.cc .
cp $SRC/src/model/layer/* .
cp $SRC/src/model/optimizer/* .
cp $SRC/src/model/loss/* .
cp $SRC/src/model/metric/* .
cp $SRC/src/model/updater/* .
cp $SRC/src/utils/* .
cp $SRC/src/model/feed_forward_net.cc .

cp -r $SRC/include/singa .
rsync -a --ignore-existing $SRC/build/include/singa/ singa/

mv cnn-bare.cc singa-cifar10.cc

rm -f cudnn*
rm -f opencl*
rm -f *.cu
rm -f *.cu
rm -f *.cl
rm -f tensor_math_opencl.h tensor_math_cuda.h
rm -f prelu.cc prelu.h rnn.cc rnn.h slice.cc slice.h dropout.cc dropout.h concat.cc concat.h \
batchnorm.cc batchnorm.h local_all_reduce.* local_updater.*
rm -rf singa/proto
