# Freebayes polish
## Update note

**2019-May-30 update**
- use Giulio's version of Freebayes which parallize by contig rather than fix bp
- get rid of bcftools consensus applet since it is now baked in Freebayes app
- add Giulio's applet for asm-stats and qv-calculator which use Arang code

**2019-May-9 update**
- explicitly add output folder in dxworkflow.json

**2019-Apr-25 update**
- prevent name collision for bcftools_consensus by setting name differently
- change default instance type for 10x longranger to mem1_ssd2_x36
