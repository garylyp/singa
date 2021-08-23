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

#include "file_reader.h"

namespace singa {

void FileReader::OpenForRead(std::string file_name) {
	if (isOpenForRead_)
		fin_.close();

	fin_.open(file_name, std::ios::in | std::ios::binary);
	CHECK(fin_.is_open()) << "Cannot read file " << file_name;

	isOpenForRead_ = true;
	isOpenForWrite_ = false;
}

void FileReader::OpenForWrite(std::string file_name) {
	if (isOpenForWrite_)
		fout_.close();

	fout_.open(file_name, std::ios::binary | std::ios::out);
	CHECK(fout_.is_open()) << "Cannot create file " << file_name;

	isOpenForRead_ = false;
	isOpenForWrite_ = true;
}

void FileReader::Close() {
	fin_.close();
	fout_.close();
}

bool FileReader::Read(std::string* key, uint8_t* bytes, size_t* size) {
	CHECK(isOpenForRead_) << "File not open for read";
	CHECK_NOTNULL(bytes);

	size_t key_len;
	fin_.read((char*)&key_len, sizeof(size_t));
	if (fin_.gcount() != sizeof(size_t))
		return false;
	char* key_str = new char[key_len];
	fin_.read(key_str, key_len);
	if (fin_.gcount() != key_len)
		return false;
	*key = std::string(key_str);
	delete key_str;
	fin_.read((char*)size, sizeof(size_t));
	if (fin_.gcount() != sizeof(size_t))
		return false;
	fin_.read((char*)bytes, *size);
	if (fin_.gcount() != *size)
		return false;

	LOG(INFO) << "Read key " << *key << " data size " << *size;

	return true;
}

std::vector<std::pair<std::string, Tensor>> FileReader::Read() {
	CHECK(isOpenForRead_) << "File not open for read";

	std::vector<std::pair<std::string, Tensor>> ret;
	fin_.seekg(0, fin_.end);
	int len = fin_.tellg();
	fin_.seekg(0, fin_.beg);
	uint8_t* buffer = new uint8_t[len];

	std::string key;
	size_t actual_size;

	while (Read(&key, buffer, &actual_size)) {
		Tensor tensor;
		tensor.FromBytes(buffer, actual_size);
		ret.push_back(std::make_pair(key, tensor));
	}

	return ret;
}

void FileReader::Write(std::string key, uint8_t* bytes, size_t size) {
	CHECK(isOpenForWrite_) << "File not open for write";

	LOG(INFO) << "Write key " << key << " data size " << size;

	size_t key_len = strlen(key.c_str()) + 1;
	fout_.write((char*)&key_len, sizeof(size_t));
	fout_.write(key.c_str(), key_len);
	fout_.write((char*)&size, sizeof(size_t));
	fout_.write((char*)bytes, size);
}

void FileReader::Write(std::string key, Tensor tensor) {
	CHECK(isOpenForWrite_) << "File not open for write";

	uint8_t* buffer = nullptr;
	size_t size;
	tensor.ToBytes(&buffer, 0, &size);
	Write(key, buffer, size);
	delete buffer;
}

} // namespace singa
