# PacBio BAM Resequencing Workflow
## Update note

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