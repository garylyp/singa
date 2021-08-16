cnn: cnn.cc
	g++ -std=c++11 -Wall -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-cuda: cnn.cc
	g++ -std=c++11 -Wall -I../../include -I../../build/include -I/usr/local/cuda/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-myfile: cnn.cc file_reader.cc
	g++ -std=c++11 -DMY_FILE_READER -Wall -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-memreader: cnn.cc mem_reader.cc myfilesnap.o
	g++ -std=c++11 -DMY_MEM_READER -Wall -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-cuda-myfile: cnn.cc file_reader.cc
	g++ -std=c++11 -DMY_FILE_READER -Wall -I../../include -I../../build/include -I/usr/local/cuda/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-bare: images.o labels.o mysnap.o mem_reader.cc cnn-bare.cc
	g++ -std=c++11 -Wall -g -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-bare-myfile: images.o labels.o myfilesnap.o mem_reader.cc cnn-bare.cc
	g++ -std=c++11 -DMY_FILE_READER -Wall -g -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-one-image: cnn-one-image.cc
	g++ -std=c++11 -DLOCAL_MAIN -Wall -g -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-one-image-filereader: cnn-one-image.cc file_reader.cc
	g++ -std=c++11 -DMY_FILE_READER -DLOCAL_MAIN -Wall -g -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-one-image-memreader: cnn-one-image.cc mem_reader.cc myfilesnap.o
	g++ -std=c++11 -DMY_MEM_READER -DLOCAL_MAIN -Wall -g -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-one-image-cuda: cnn-one-image.cc
	g++ -std=c++11 -DLOCAL_MAIN -Wall -g -I../../include -I../../build/include -I/usr/local/cuda/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

cnn-hadoop-pipes: cnn-hadoop-pipes.cc cnn-one-image.cc
	g++ -std=c++11 -Wall -g -I../../include -I../../build/include -I/usr/local/cuda/include -I$(WORKING_DIR)/hadoop-2.10.1/include $^ -o $@ -L../../build/lib -L/usr/local/lib -L$(WORKING_DIR)/hadoop-2.10.1/lib/native -L$(JAVA_HOME)/jre/lib/amd64/server -lsinga -lprotobuf -lglog -lhdfs -lhadooppipes -lpthread -ljvm -lcrypto

vgg: vgg.cc
	g++ -std=c++11 -Wall -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

vgg-cuda: vgg.cc
	g++ -std=c++11 -Wall -I../../include -I../../build/include -I/usr/local/cuda/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

split-dataset-images: split-dataset-images.c
	gcc -Wall $^ -o $@

objects:  
	objcopy -I binary -O elf64-x86-64 -B i386:x86-64 images.bin images.o
	objcopy -I binary -O elf64-x86-64 -B i386:x86-64 mysnap.bin mysnap.o
	objcopy -I binary -O elf64-x86-64 -B i386:x86-64 myfilesnap.bin myfilesnap.o
	
cnn-one-image-bare: cnn-one-image-bare.cc mem_reader.cc images.o myfilesnap.o
	g++ -g -I../../build/include -I../../include $^ -o $@ -L../../build/lib -lsinga -lprotobuf

convert: convert.cc file_reader.cc
	g++ -DMY_FILE_READER -I../../build/include -I../../include $^ -o $@ -L../../build/lib -lsinga -lprotobuf

lenet: lenet.cc
	g++ -std=c++11 -Wall -I../../include -I../../build/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

lenet-cuda: lenet.cc
	g++ -std=c++11 -Wall -I../../include -I../../build/include -I/usr/local/cuda/include $^ -o $@ -L../../build/lib -L/usr/local/lib -lsinga -lprotobuf -lglog

clean:
	rm -f cnn cnn-cuda
