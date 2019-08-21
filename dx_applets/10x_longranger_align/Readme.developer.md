Arang Rhie 4:37 PM
@chai_temp are you using override.json style fix for longranger?

Arang Rhie 4:37 PM
{
"ALIGNER_CS.ALIGNER._LINKED_READS_ALIGNER.BARCODE_AWARE_ALIGNER": { "chunk.mem_gb": 48 },
"ALIGNER_CS.ALIGNER._LINKED_READS_ALIGNER.MERGE_POS_BAM": { "join.mem_gb": 48 },
"ALIGNER_CS.ALIGNER._REPORTER.FILTER_BARCODES": { "join.mem_gb": 48 },
"ALIGNER_CS.ALIGNER._REPORTER.REPORT_LENGTH_MASS": { "chunk.mem_gb": 32 }
}

So if the instance if out of memory, the user not only have to increase the instance type, but also have tell the long ranger to take up more memory as well. Is that correct?
