#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 file partSizeInMb";
    exit 0;
fi

file=$1

if [ ! -f "$file" ]; then
    echo "Error: $file not found." 
    exit 1;
fi

partSizeInMb=8  # Default: 8 Mb. Change this if you used different part sizes.
fileSizeInMb=$(du -b "$file" | cut -f 1 | awk '{printf("%.0f\n", $1/1048576)}')
parts=$((fileSizeInMb / partSizeInMb))
if [[ $((fileSizeInMb % partSizeInMb)) -gt 0 ]]; then
    parts=$((parts + 1));
fi

if [ $parts -eq 0 ]; then
	md5sum $file | awk '{print $1}'
	exit
fi

checksumFile=$(mktemp -t s3md5.XXXXXX)

for (( part=0; part<$parts; part++ ))
do
    skip=$((partSizeInMb * part))
    partSizeInB=$((partSizeInMb * 1000000))
    $(dd bs=1024k count=$partSizeInMb skip=$skip if="$file" 2>/dev/null | md5sum >>$checksumFile)
done

echo "$(xxd -r -p $checksumFile | md5sum | awk '{print $1}')-$parts"
rm $checksumFile
