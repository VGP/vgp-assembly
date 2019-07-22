<!-- dx-header -->
# proc10xG: trimmer of GEM barcodes - VGP

This apps removes the GEM barcodes from paired 10x reads. It accepts arrays of files.
Filenames need to be in the form <Specied_ID>_S1_L001_R1_001.fastq.gz and paired.

You can run it on small-to-large instances, depending on the number of files to process.
It will use a core for each pair in parallel, so if you have 8 paired files, you should use an instance with 4 cores.  Don't forget to increase disk space accordingly

<!-- /dx-header -->

<!-- proc10xG: trimmer of GEM barcodes - VGP

This apps removes the GEM barcodes from paired 10x reads. It accepts arrays of files.
Filenames need to be in the form <Specied_ID>_S1_L001_R1_001.fastq.gz and paired.

You can run it on small-to-large instances, depending on the number of files to process.
It will use a core for each pair in parallel, so if you have 8 paired files, you should use an instance with 4 cores.  Don't forget to increase disk space accordingly

-->
