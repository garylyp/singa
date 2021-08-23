<!--
    Licensed to the Apache Software Foundation (ASF) under one
    or more contributor license agreements.  See the NOTICE file
    distributed with < this work for additional information
    regarding copyright ownership.  The ASF licenses this file
    to you under the Apache License, Version 2.0 (the
    "License"); you may not use this file except in compliance
    with the License.  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an
    "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
    KIND, either express or implied.  See the License for the
    specific language governing permissions and limitations
    under the License.
-->

![Logo](doc/_static/singa.png)

# Apache SINGA

![Native Ubuntu build status](https://github.com/apache/singa/workflows/Native-Ubuntu/badge.svg)
![Native Mac build status](https://github.com/apache/singa/workflows/Native-MacOS/badge.svg)
![conda build status](https://github.com/apache/singa/workflows/conda/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/apache-singa/badge/?version=latest)](https://apache-singa.readthedocs.io/en/latest/?badge=latest)
![License](http://img.shields.io/:license-Apache%202.0-blue.svg)
[![Follow Apache SINGA on Twitter](https://img.shields.io/twitter/follow/apachesinga.svg?style=social&label=Follow)](https://twitter.com/ApacheSinga)
[![Docker pulls](https://img.shields.io/docker/pulls/apache/singa.svg)](https://hub.docker.com/r/apache/singa/)

Distributed deep learning system

[http://singa.apache.org](http://singa.apache.org)

## Quick Start

* [Installation](http://singa.apache.org/docs/installation/)
* [Examples](examples)

## Issues

* [JIRA tickets](https://issues.apache.org/jira/browse/SINGA)

## Code Analysis:

![LGTM C++ Grade](https://img.shields.io/lgtm/grade/cpp/github/apache/singa)
![LGTM Python Grade](https://img.shields.io/lgtm/grade/python/github/apache/singa)
[![codecov](https://codecov.io/gh/apache/singa/branch/master/graph/badge.svg)](https://codecov.io/gh/apache/singa)

[![Stargazers over time](https://starchart.cc/apache/singa.svg)](https://starchart.cc/apache/singa)

## Mailing Lists

* [Development Mailing List](mailto:dev-subscribe@singa.incubator.apache.org) ([Archive](http://mail-archives.apache.org/mod_mbox/singa-dev/))
* [Commits Mailing List](mailto:commits-subscribe@singa.incubator.apache.org) ([Archive](http://mail-archives.apache.org/mod_mbox/singa-commits/))

## This Branch

This branch contains the code to build fat executables to be run on a posit hardware.

How to run the code on an x86/64 host:

```
cd tool/posit
./copy-srcs.sh
cd ../../build_posit
make -j4
```

You should train the cnn cifar-10 model and obtain a snapshot:

```
cd examples/cifar10
make
./download_data.py bin
./cnn
```

This should create two file ``mysnap.desc`` and ``mysnap.bin``.

```
cd build_posit
cp ../examples/cifar10/mysnap.desc .
cp ../examples/cifar10/mysnap.bin .
cp -r cp ../examples/cifar10/cifar-10-batches-bin .
./cnn
```

The result should be:

```
...
Accuracy: 0.4774
```

Build docker image with posit support:

```
docker build --rm=false /home/dumi/git/incubator-singa/tool/docker/devel/ubuntu/cuda10-posit -t dumi/singa:devel-cuda10-cudnn7-posit
```

TBD: how to compile and run on the posit hardware.
