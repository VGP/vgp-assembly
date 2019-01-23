# FALCON and FALCON-Unzip Assembly
## Update note

**2018-Jan-23 update**
- update to Falcon falcon-2018.31.08-03.06
- change default instance type Unzip stage 3 to mem4_ssd1_x128
- fix jellyfish/genomescope incorrect link from Daligner stage2 to stage1
- make dazzler name different among standard, tanmask, and repmask
- upgrade Unzip stage 5 to support V6 chemistry using SMRTlink 6.0.0.47841
- change parameter for Falcon stage 1 from -k16 -e0.70 -s1000 -l1000 -h64 -w7 to -k14 -e0.75 -s100 -l2500 -h240 -w8 and for stage 2 from -k20 -e.96 -s1000 -t32 -l2500 -h256 to -k24 -e.90 -s100 -l1000 -h600. 

Backward incompatibility note:
- dazzler db created by this version of workflow with _tan and _mask cannot be used with older version of "HPC REPmask", "Daligner for PacBio's Falcon", or "Falcon_asm" apps/applets. However, output from old version of apps/applets would work fine for the new version (this version).