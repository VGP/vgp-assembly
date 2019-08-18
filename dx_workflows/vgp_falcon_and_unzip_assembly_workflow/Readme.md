# FALCON and FALCON-Unzip Assembly
## Update note

**2019-Aug-18 update**
- debug unzip_4_haplotype_polish
- increase timeout for unzip_3_haplotype_track_reads due to recent long running of mPanTro1

**2019-Jul-31 update**
- update bam_to_fasta and unzip_4_haplotype_polish to SMRTlink7 to support Sequel II
- change tan and rep mask second stage to newer version which we forgot in 2019-Feb-15

**2019-Feb-26 update**
- change default instance type for Align and Phase

**2019-Feb-15 update**
- start using TANmask2.0.0 and REPmask2.0.1 after testing on 3 large genomes
- fix discrepancy of unzip_stage_0 app and github. The old version miss logic to handle missing ref/read. New version is unzip_stage_0 track_read 1.1.0
- add calculate read length distribution for primary and haplotig

**2019-Jan-23patch1 update**
- new TANmask and REPmask is very buggy, so I roll back TAN and REPmask from 2.0.0 to 1.9.3. That version has difference name for standard, tanmask, and repmask, but use old Falcon code.

**2019-Jan-23 update**
- update to Falcon falcon-2018.31.08-03.06
- change default instance type Unzip stage 3 to mem4_ssd1_x128
- fix jellyfish/genomescope incorrect link from Daligner stage2 to stage1
- make dazzler name different among standard, tanmask, and repmask
- upgrade Unzip stage 5 to support V6 chemistry using SMRTlink 6.0.0.47841
- change parameter for Falcon stage 1 from -k16 -e0.70 -s1000 -l1000 -h64 -w7 to -k14 -e0.75 -s100 -l2500 -h240 -w8 and for stage 2 from -k20 -e.96 -s1000 -t32 -l2500 -h256 to -k24 -e.90 -s100 -l1000 -h600. 

Backward incompatibility note:
- dazzler db created by this version of workflow with _tan and _mask cannot be used with older version of "HPC REPmask", "Daligner for PacBio's Falcon", or "Falcon_asm" apps/applets. However, output from old version of apps/applets would work fine for the new version (this version).
