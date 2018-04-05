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

| species_name	| species_id	| tech_count	| pacbio_subreads	| pacbio_scrubs	| 10x	| bionano_tgh	| bionano_dls	| bionano_bnx	| bionano_cmap	| hic |
| :---------- | :---------- | :---------- | :---------- | :---------- | :----- | :----- | :----- | :----- | :----- | :----- |
| Calypte_anna	| bCalAnn1	| 4	| 189	| 0	| 1	| O	| O	| O	| O	| arima |
| Rhinolophus_ferrumequinum	| mRhiFer1	| 3	| 25	| 25	| 12	| O	| X	| O	| O	| phase |
| Anabas_testudineus	| fAnaTes1	| 3	| 8	| 0	| 4	| O	| X	| O	| O	| arima |
| Mastacembelus_armatus	| fMasArm1	| 3	| 8	| 0	| 4	| O	| X	| O	| O	| arima |
| Sparus_aurata	| fSpaAur1	| 2	| 10	| 0	| 12	| X	| X	| X	| X	|  |
| Rhinatrema_bivittatum	| aRhiBiv1	| 2	| 54	| 0	| 12	| O	| X	| O	| O	|  |
| Gouania_willdenowi	| fGouWil2	| 2	| 16	| 0	| 4	| O	| X	| O	| O	|  |
| Aquila_chrysaetos	| bAquChr1	| 2	| 16	| 0	| 12	| X	| X	| X	| X	|  |
| Denticeps_clupeoides	| fDenClu1	| 2	| 8	| 0	| 4	| X	| X	| X	| X	|  |
| Astatotilapia_calliptera	| fAstCal1	| 2	| 67	| 0	| 4	| X	| X	| O	| O	|  |
| Streptopelia_turtur	| bStrTur1	| 2	| 18	| 0	| 12	| X	| X	| X	| X	|  |
| Erpetoichthys_calabaricus	| fErpCal1	| 2	| 53	| 0	| 8	| X	| X	| X	| X	|  |
| Archocentrus_centrarchus	| fArcCen1	| 2	| 18	| 0	| 0	| X	| X	| X	| X	| phase |
| Zeus_faber	| fZeuFab1	| 2	| 27	| 0	| 4	| X	| X	| X	| X	|  |
| Cottoperca_gobio	| fCotGob3	| 2	| 10	| 0	| 4	| O	| X	| O	| O	|  |
| Gopherus_evgoodei	| rGopEvg1	| 1	| 19	| 19	| 0	| X	| X	| X	| X	|  |
| Amblyraja_radiata	| sAmbRad1	| 1	| 43	| 42	| 0	| O	| X	| O	| X	|  |
