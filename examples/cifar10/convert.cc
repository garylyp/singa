#include "singa/io/snapshot.h"
#include "file_reader.h"

using namespace singa;

int main() {
	Snapshot snap("mysnap", Snapshot::kRead, 100);
	vector<std::pair<std::string, Tensor>> params = snap.Read();
	FileReader fsnap;
	fsnap.OpenForWrite("myfilesnap.bin");
	
	for (auto it = params.end() - 1; ; it--) {
		fsnap.Write((*it).first, (*it).second);
		if (it == params.begin())
			break;
	}
	fsnap.Close();
	return 0;
}
