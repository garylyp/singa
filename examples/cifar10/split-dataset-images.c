#include <stdio.h>

#define BIN_FILE "cifar-10-batches-bin/test_batch.bin"
#define OUT_PREFIX "cifar-10-split-bin"
#define OUT_LABELS "cifar-10-split-bin/labels.txt"

#define N_IMG	10000
#define IMG_DIM	32
#define IMG_BIN_SIZE IMG_DIM * IMG_DIM * 3 + 1

int main(int argc, char** argv) {
	char fn_buff[64];
	char img_buff[IMG_BIN_SIZE];
	FILE* fin = fopen(BIN_FILE, "rb");
	if (!fin)
		return -1;
	FILE* flbl = fopen(OUT_LABELS, "wt");
	if (!flbl) {
		fclose(fin);
		return -1;
	}
	for (int i = 0; i < N_IMG; i++) {
		int nb = fread(img_buff, sizeof(char), IMG_BIN_SIZE, fin);
		if (nb != IMG_BIN_SIZE) {
			printf("Error reading %d bytes index %d!", IMG_BIN_SIZE, i);
			goto finish;
		}
		sprintf(fn_buff, "%s/%04d.bin", OUT_PREFIX, i);
		FILE* fout = fopen(fn_buff, "wb");
		if (!fout) {
			printf("Error creating file %s\n", fn_buff);
			goto finish;
		}
		nb = fwrite(img_buff, sizeof(char), IMG_BIN_SIZE, fout);
		if (nb != IMG_BIN_SIZE) {
			printf("Error writing %d bytes index %d!", IMG_BIN_SIZE, i);
			fclose(fout);
			goto finish;
		}
		fclose(fout);
		fprintf(flbl, "%04d.bin %d\n", i, (int)img_buff[0]);
	}
	finish:
	fclose(flbl);
	fclose(fin);
	return 0;
}
