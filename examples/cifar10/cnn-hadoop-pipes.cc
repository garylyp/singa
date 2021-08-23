/*
 * SingaCNN Hadoop Pipes C++.
 *
 * Based on Apache Singa 2.0
 *
 * Copyright (C) June 2020  Dumi Loghin
 */

#include <sys/time.h>
#include <iostream>

#include "cnn-one-image.h"

#include "hdfs.h"
#include "Pipes.hh"
#include "TemplateFactory.hh"
#include "StringUtils.hh"

using namespace std;

#define StrDatasetKey 			"hdfs.cifar10.datasetpath"
#define SingleReduceKey			"image_label_prediction"
#define ReduceOutputPredictions	"predictions"
#define ReduceOutputAccuracy	"accuracy"

#define DIM 3073				// 1 + 32 * 32 * 3

class SingaCNNMapper : public HadoopPipes::Mapper {

private:
	string datasetpath;
	hdfsFS dfs;

public:
	SingaCNNMapper(HadoopPipes::TaskContext& context) {
		singa::InitChannel(nullptr);
		HadoopPipes::JobConf* conf = (HadoopPipes::JobConf*)context.getJobConf();
		if (conf != NULL) {
			datasetpath = conf->get(StrDatasetKey);
		}
		else {
			cerr << "Fatal: no JobConf found!" << endl;
			exit(-1);
		}
		dfs = hdfsConnectAsUser("hdfs://localhost", 54310, "hadoop");
		if (!dfs) {
			cerr << "Error connecting to HDFS!" << endl;
			exit(-1);
		}
	}

	void map(HadoopPipes::MapContext& context) {
		string val = context.getInputValue();
		int pos = val.find(' ');
		string key = val.substr(0, pos);
		string path = datasetpath + "/" + key;
		hdfsFile img_file = hdfsOpenFile(dfs, path.c_str(), O_RDONLY, 0, 0, 0);
		if(!img_file) {
			cerr << "Failed to open file for reading " << path << endl;
			return;
		}
		char buff[DIM];
		tSize nb = hdfsRead(dfs, img_file, (void*)buff, DIM);
		if (nb != DIM) {
			cerr << "Error reading image: " << nb << " bytes read, " << DIM << " expected!" << endl;
			return;
		}
		hdfsCloseFile(dfs, img_file);
		int idx = singa::EvalFromBuffer(buff);
		context.emit(SingleReduceKey, val.append(" ").append(HadoopUtils::toString(idx)));
	}
};

class SingaCNNReducer : public HadoopPipes::Reducer {
public:
	SingaCNNReducer(HadoopPipes::TaskContext& context) {}

	void reduce(HadoopPipes::ReduceContext& context) {
		string key = context.getInputKey();

		while (context.nextValue()) {
			string key = context.getInputKey();
			int n = 0;
			int match = 0;

			// iterate through all the points
			vector<long> labels;
			while (context.nextValue()) {
				string val = context.getInputValue();
				labels.clear();
				int pos1 = val.find(' ');
				int pos2 = val.find(' ', pos1+1);
				long label = strtol(val.substr(pos1+1, pos2).c_str(), NULL, 10);
				long prediction = strtol(val.substr(pos2+1).c_str(), NULL, 10);
				if (label == prediction)
					match++;
				n++;
				context.emit(ReduceOutputPredictions, val);
			}
			char buff[8];
			sprintf(buff, "%4.2f%%", 100.0 * match / n);
			context.emit(ReduceOutputAccuracy, buff);
		}
	}
};

int main(int argc, char** argv) {
	//start a map/reduce job
	return HadoopPipes::runTask(HadoopPipes::TemplateFactory<SingaCNNMapper, SingaCNNReducer>());
}
