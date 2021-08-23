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

#ifndef SINGA_CIFAR10_OBJECTS_H_
#define SINGA_CIFAR10_OBJECTS_H_

#define USE_FLOAT_IMAGE

#ifdef USE_FLOAT_IMAGE
extern unsigned char _binary_images_bin_start;
#else
extern unsigned char _binary_images_int_bin_start;
#endif

extern unsigned char _binary_labels_bin_start;
extern unsigned char _binary_mysnap_bin_start;
extern unsigned char _binary_mysnap_bin_end;
extern unsigned char _binary_mysnap_bin_size;

extern unsigned char _binary_myfilesnap_bin_start;
extern unsigned char _binary_myfilesnap_bin_end;
extern unsigned char _binary_myfilesnap_bin_size;

#endif
