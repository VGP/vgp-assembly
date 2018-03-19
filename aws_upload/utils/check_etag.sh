#!/bin/bash


######################################################################
#
#  This software is based on:
#    'A Bash script to compute ETag values for S3 multipart uploads on OS X.'
#    (https://gist.github.com/emersonf/7413337)
#
#  Modifications are a 'United States Government Work', and
#  are released in the public domain.
#  
#
#  PUBLIC DOMAIN NOTICE
#
#  This software is "United States Government Work" under the terms of the United
#  States Copyright Act. It was written as part of the authors' official duties
#  for the United States Government and thus cannot be copyrighted. This software
#  is freely available to the public for use without a copyright
#  notice. Restrictions cannot be placed on its present or future use.
#
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and associated data, the National Human Genome
#  Research Institute (NHGRI), National Institutes of Health (NIH) and the
#  U.S. Government do not and cannot warrant the performance or results that may
#  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
#  disclaim all warranties as to performance, merchantability or fitness for any
#  particular purpose.
#
#  Please cite the authors in any work or product based on this material.
######################################################################

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
