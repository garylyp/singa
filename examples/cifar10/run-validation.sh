#!/bin/bash

cd cifar-10-split-bin
FILES=`ls *.bin`
cd ..

OFILE="validation-labels.txt"
cat /dev/null > $OFILE

for FILE in $FILES; do
	cp cifar-10-split-bin/$FILE images.bin
	make objects > /dev/null 2>&1
	rm cnn-one-image-bare
	make cnn-one-image-bare > /dev/null 2>&1
	./cnn-one-image-bare > cifar-10-split-bin/$FILE.txt 2>&1
	LBL=`cat cifar-10-split-bin/$FILE.txt | grep Label | cut -d ' ' -f 7`
	echo "$FILE $LBL" >> $OFILE
done
echo "Done."