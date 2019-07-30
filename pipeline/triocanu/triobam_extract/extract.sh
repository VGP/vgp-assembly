#!/bin/bash

BAM_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bam.uniq)

BAM_FILE=$(grep ${BAM_NAME} bam.ls)

if ! [[ -e "filtered/${BAM_NAME}.subreads.bam" ]]; then

	java -jar ${PICARD_PATH}/picard.jar FilterSamReads I=${BAM_FILE} O=filtered/${BAM_NAME}.subreads.bam READ_LIST_FILE=reads/${BAM_NAME}.ls VALIDATION_STRINGENCY=LENIENT FILTER=includeReadList
	
fi