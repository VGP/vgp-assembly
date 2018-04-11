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

### Status as of 2018-04-10

O: Available

X: Not Available

Total num. objects under s3://genomeark/species: 2,463


| species_name	| species_id	| tech_count	| pacbio_subreads	| pacbio_scrubs	| 10x	| bionano_tgh	| bionano_dls	| bionano_bnx	| bionano_cmap	| hic |
| :---------- | :---------- | :---------- | :---------- | :---------- | :----- | :----- | :----- | :----- | :----- | :----- |
| Calypte_anna	| bCalAnn1	| 4	| 189	| 0	| 1	| O	| O	| O	| O	| arima |
| Anabas_testudineus	| fAnaTes1	| 4	| 8	| 8	| 4	| O	| X	| O	| O	| arima |
| Rhinolophus_ferrumequinum	| mRhiFer1	| 4	| 25	| 25	| 12	| O	| X	| O	| O	| phase |
| Mastacembelus_armatus	| fMasArm1	| 4	| 8	| 0	| 4	| O	| X	| O	| O	| arima |
| Phyllostomus_discolor	| mPhyDis1	| 3	| 43	| 43	| 0	| O	| X	| O	| X	| arima |
| Rhinatrema_bivittatum	| aRhiBiv1	| 3	| 54	| 0	| 12	| O	| X	| O	| O	|  |
| Taeniopygia_guttata	| bTaeGut2	| 3	| 17	| 17	| 0	| O	| X	| O	| X	| arima |
| Turdus_merula	| bTurMer1	| 2	| 11	| 11	| 0	| O	| O	| O	| X	|  |
| Gouania_willdenowi	| fGouWil2	| 3	| 16	| 16	| 4	| O	| X	| O	| O	|  |
| Amblyraja_radiata	| sAmbRad1	| 3	| 50	| 50	| 0	| O	| X	| O	| X	| arima |
| Lynx_canadensis	| mLynCan4	| 3	| 35	| 35	| 0	| O	| X	| O	| X	| arima |
| Cottoperca_gobio	| fCotGob3	| 3	| 10	| 10	| 4	| O	| X	| O	| O	|  |
| Dendrocygna_viduata	| bDenVid1	| 2	| 10	| 10	| 0	| O	| X	| O	| X	|  |
| Sparus_aurata	| fSpaAur1	| 2	| 10	| 0	| 12	| X	| X	| X	| X	|  |
| Spatula_cyanoptera	| bSpaCya1	| 2	| 18	| 18	| 0	| O	| X	| O	| X	|  |
| Aquila_chrysaetos	| bAquChr1	| 2	| 16	| 0	| 12	| X	| X	| X	| X	|  |
| Denticeps_clupeoides	| fDenClu1	| 2	| 8	| 8	| 4	| X	| X	| X	| X	|  |
| Astatotilapia_calliptera	| fAstCal1	| 2	| 67	| 0	| 4	| X	| X	| O	| O	|  |
| Branta_leucopsis	| bBraLeu1	| 2	| 16	| 16	| 0	| O	| X	| O	| X	|  |
| Archocentrus_centrarchus	| fArcCen1	| 2	| 18	| 0	| 0	| X	| X	| X	| X	| phase |
| Zeus_faber	| fZeuFab1	| 2	| 27	| 0	| 4	| X	| X	| X	| X	|  |
| Gopherus_evgoodei	| rGopEvg1	| 2	| 19	| 19	| 0	| X	| O	| O	| X	|  |
| Strigops_habroptilus	| bStrHab1	| 2	| 0	| 0	| 0	| X	| O	| O	| X	| arima |
| Streptopelia_turtur	| bStrTur1	| 2	| 18	| 0	| 12	| X	| X	| X	| X	|  |
| Erpetoichthys_calabaricus	| fErpCal1	| 2	| 53	| 53	| 8	| X	| X	| X	| X	|  |
| Geothlypis_trichas	| bGeoTri1	| 2	| 0	| 0	| 0	| O	| X	| O	| X	| arima |
| Taeniopygia_guttata	| bTaeGut1	| 1	| 0	| 0	| 0	| X	| X	| X	| X	| arima |
| Nomonyx_dominicus	| bNomDom1	| 1	| 0	| 0	| 0	| O	| X	| O	| X	|  |
