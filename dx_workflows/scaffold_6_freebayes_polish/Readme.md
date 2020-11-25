# Freebayes polish
## Update note

**2020-Nov-25 update**
- use only QV function for merqury due to unknown bug for using full calculation

**2020-Oct-21 update**
- update merqury Docker 1.1 which fix trio calculation and give full report

**2020-Apr-21 update**
- update to merqury Docker 0.0.2 (the Apr 18 was mistakenly use old docker)

**2020-Apr-18 update**
- add merqury QV 0.0.2

**2019-Oct-07 update**
- use mem3_ssd1_v2_x48 for 10x align instance

**2019-Sep-03 update**
- use fastq.gz directly for longranger 10x align

**2019-Jul-08 update**
- update QV to Giulio new version

**2019-Jun-27 update**
- add extension guide for asm-qv

**2019-Jun-26 update**
- add extension guide for freebayes

**2019-Jun-23 update**
- fix freebayes quote issue by hardcode inclusion statement

**2019-Jun-22 update**
- change qv applet to use larger storage and longer timeout

**2019-Jun-21 update**
- added and disable distributed run in Freebayes
- use summary.csv to guide skip-coverage 

**2019-May-30 update**
- use Giulio's version of Freebayes which parallize by contig rather than fix bp
- get rid of bcftools consensus applet since it is now baked in Freebayes app
- add Giulio's applet for asm-stats and qv-calculator which use Arang code

**2019-May-9 update**
- explicitly add output folder in dxworkflow.json

**2019-Apr-25 update**
- prevent name collision for bcftools_consensus by setting name differently
- change default instance type for 10x longranger to mem1_ssd2_x36
