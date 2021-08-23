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

#ifndef SINGA_MEM_READER_H_
#define SINGA_MEM_READER_H_

#include <string>
#include <cstring>
#include <vector>

#include "singa/core/tensor.h"

namespace singa {

class MemReader {
 public:

  MemReader(char* src, int size);

  void SeekToFirst();

  bool Read(std::string* key, std::string* value);

  bool Read(std::string* key, char** bytes, size_t* size);

  std::vector<std::pair<std::string, Tensor>> Read();

 protected:
  bool ReadField(std::string* content);

 private:
  /// internal buffer
  char* buf_;
  /// offset inside the buf_
  int offset_;
  /// bytes in buf_
  int bufsize_;
  /// magic word
  const char kMagicWord[2] = {'s', 'g'};
};

} // end namespace singa

#endif // SINGA_MEM_READER_H_
