# Freebayes polish
## Update note

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
