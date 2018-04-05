# Data Sharing Status of GenomeArk

This script retrieves objects from the public genomeark and returns a short summary.

### Compile
```
javac -cp .:third-party/lib/* UpdateGenomeArkStatus.java
```

### Run
```
java -cp .:third-party/lib/* UpdateGenomeArkStatus
```

### Output format option
To create a .md compatable table, run with the option md.
By default, data will be printed in tab-delimited format.
```
java -cp .:third-party/lib/* UpdateGenomeArkStatus -md
```

### Status as of 2018-04-05

| species_name	| species_id	| pacbio_subreads	| pacbio_scrubs	| 10x	| bionano_tgh	| bionano_dls	| bionano_bnx	| bionano_cmap	| hic |
| :---------- | :---------- | :---------- | :---------- | :----- | :----- | :----- | :----- | :----- | :----- |
| Sparus_aurata	| fSpaAur1	| 20	| 0	| 12	| X	| X	| X	| X	|  |
| Gopherus_evgoodei	| rGopEvg1	| 38	| 38	| 0	| X	| X	| X	| X	|  |
| Rhinatrema_bivittatum	| aRhiBiv1	| 108	| 0	| 12	| O	| X	| O	| O	|  |
| Rhinolophus_ferrumequinum	| mRhiFer1	| 50	| 25	| 12	| O	| X	| O	| O	| phase |
| Gouania_willdenowi	| fGouWil2	| 32	| 0	| 4	| O	| X	| O	| O	|  |
| Aquila_chrysaetos	| bAquChr1	| 32	| 0	| 12	| X	| X	| X	| X	|  |
| Denticeps_clupeoides	| fDenClu1	| 16	| 0	| 4	| X	| X	| X	| X	|  |
| Astatotilapia_calliptera	| fAstCal1	| 67	| 0	| 4	| X	| X	| O	| O	|  |
| Amblyraja_radiata	| sAmbRad1	| 36	| 34	| 0	| O	| X	| O	| X	|  |
| Calypte_anna	| bCalAnn1	| 378	| 0	| 1	| O	| O	| O	| O	| arima |
| Streptopelia_turtur	| bStrTur1	| 36	| 0	| 12	| X	| X	| X	| X	|  |
| Erpetoichthys_calabaricus	| fErpCal1	| 106	| 0	| 8	| X	| X	| X	| X	|  |
| Archocentrus_centrarchus	| fArcCen1	| 18	| 0	| 0	| X	| X	| X	| X	| phase |
| Zeus_faber	| fZeuFab1	| 54	| 0	| 4	| X	| X	| X	| X	|  |
| Cottoperca_gobio	| fCotGob3	| 20	| 0	| 4	| O	| X	| O	| O	|  |
| Anabas_testudineus	| fAnaTes1	| 16	| 0	| 4	| O	| X	| O	| O	| arima |
| Mastacembelus_armatus	| fMasArm1	| 16	| 0	| 4	| O	| X	| O	| O	| arima |
