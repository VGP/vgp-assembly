# PacBio BAM Resequencing Workflow
## Update note

**2020-Apr-21 update**
- update to merqury Docker 0.0.2 (the Apr 18 was mistakenly use old docker)

**2020-Apr-18 update**
- add QV measurement 0.0.2

**2020-Apr-09 update**
- update to app run_polish and app_polish version pbgcpp1.9.0.1 (internal version correspond to pbgcpp1.9.0) which has proper parameter passing to Docker image and have default for minCov10 maxCov120 confidence20

**2020-Apr-01 update**
- update to pbgcpp1.9.0

**2019-Sep-14 update**
- change default instance type for polishing

**2019-Sep-12 update**
- remove invalid instance type option

**2019-Sep-11 update**
- adjust instance type

**2019-Sep-10 update**
- use V2
- update to SMRTlink7.0.1

**2019-Apr-25 update**
- use bam_slicing and use default instance for slicer_instance_type=mem1_hdd2_x32 and polish_instance_tpye=mem1_ssd1_x32

**2019-Apr-23 update**
- change default instance type for polishing to mem3_ssd1_x32 to reduce out of storage error
- get rid of mapping report applet from the workflow since it fails frequently and no on use it.

**2019-Apr-18 update**
- change mapper from blasr to Minimap2
- change version of SMRTlink polish to 6.0
- get rid of resequencing report applet from the workflow since it is broken.