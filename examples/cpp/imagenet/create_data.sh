#!/usr/bin/env sh
#
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 

if ! [ -e "createdata" ]; then
	g++ -O3 -I../../../include -I../../../build/include ilsvrc12.cc -o createdata -L../../../build/lib -lsinga -lprotobuf -lpthread -lopencv_core -lopencv_imgcodecs -lopencv_imgproc
fi

export LD_LIBRARY_PATH=../../../build/lib

./createdata -trainlist "imagenet/label/train.txt" -trainfolder "imagenet/ILSVRC2012_img_train" \
  -testlist "imagenet/label/val.txt" -testfolder "imagenet/ILSVRC2012_img_val" -outdata "imagenet_data" -filesize 1280
