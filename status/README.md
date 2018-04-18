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

### Status as of 2018-04-18

O: Available

X: Not Available

Total num. objects under s3://genomeark/species: 9,282

| species_name	| species_id	| tech_count	| pacbio_subreads	| pacbio_scrubs	| 10x	| bionano_tgh	| bionano_dls	| bionano_bnx	| bionano_cmap	| hic |
| :---------- | :---------- | :---------- | :---------- | :---------- | :----- | :----- | :----- | :----- | :----- | :----- |
| Calypte_anna	| bCalAnn1	| 4	| 189	| 0	| 1	| O	| O	| O	| O	| arima |
| Phyllostomus_discolor	| mPhyDis1	| 4	| 43	| 43	| 8	| O	| X	| O	| X	| arima |
| Taeniopygia_guttata	| bTaeGut2	| 4	| 20	| 20	| 8	| O	| X	| O	| X	| arima |
| Taeniopygia_guttata	| bTaeGut1	| 4	| 118	| 0	| 8	| O	| X	| O	| X	| arima |
| Anabas_testudineus	| fAnaTes1	| 4	| 8	| 8	| 4	| O	| X	| O	| O	| arima |
| Rhinolophus_ferrumequinum	| mRhiFer1	| 4	| 25	| 25	| 12	| O	| X	| O	| O	| phase |
| Gouania_willdenowi	| fGouWil2	| 4	| 16	| 16	| 4	| O	| X	| O	| O	| arima |
| Amblyraja_radiata	| sAmbRad1	| 4	| 50	| 50	| 8	| O	| X	| O	| X	| arima |
| Lynx_canadensis	| mLynCan4	| 4	| 35	| 35	| 8	| O	| X	| O	| X	| arima |
| Mastacembelus_armatus	| fMasArm1	| 4	| 8	| 8	| 4	| O	| X	| O	| O	| arima |
| Sparus_aurata	| fSpaAur1	| 3	| 10	| 10	| 12	| X	| X	| X	| X	| arima |
| Rhinatrema_bivittatum	| aRhiBiv1	| 3	| 54	| 54	| 12	| O	| X	| O	| O	|  |
| Denticeps_clupeoides	| fDenClu1	| 3	| 8	| 8	| 4	| X	| X	| X	| X	| arima |
| Turdus_merula	| bTurMer1	| 2	| 11	| 11	| 0	| O	| O	| O	| X	|  |
| Strigops_habroptilus	| bStrHab1	| 3	| 126	| 0	| 0	| X	| O	| O	| X	| arima |
| Streptopelia_turtur	| bStrTur1	| 3	| 18	| 18	| 12	| X	| X	| X	| X	| arima |
| Cottoperca_gobio	| fCotGob3	| 3	| 10	| 10	| 4	| O	| X	| O	| O	|  |
| Geothlypis_trichas	| bGeoTri1	| 3	| 0	| 0	| 8	| O	| X	| O	| X	| arima |
| Dendrocygna_viduata	| bDenVid1	| 2	| 10	| 10	| 0	| O	| X	| O	| X	|  |
| Spatula_cyanoptera	| bSpaCya1	| 2	| 18	| 18	| 0	| O	| X	| O	| X	|  |
| Aquila_chrysaetos	| bAquChr1	| 2	| 16	| 16	| 12	| X	| X	| X	| X	|  |
| Astatotilapia_calliptera	| fAstCal1	| 2	| 67	| 0	| 4	| X	| X	| O	| O	|  |
| Branta_leucopsis	| bBraLeu1	| 2	| 16	| 16	| 0	| O	| X	| O	| X	|  |
| Archocentrus_centrarchus	| fArcCen1	| 2	| 18	| 0	| 0	| X	| X	| X	| X	| phase |
| Zeus_faber	| fZeuFab1	| 2	| 27	| 27	| 4	| X	| X	| X	| X	|  |
| Takifugu_rubripes	| fTakRub1	| 2	| 4	| 4	| 0	| X	| X	| X	| X	| arima |
| Gopherus_evgoodei	| rGopEvg1	| 2	| 19	| 19	| 0	| X	| O	| O	| X	|  |
| Erpetoichthys_calabaricus	| fErpCal1	| 2	| 53	| 53	| 8	| X	| X	| X	| X	|  |
| Ornithorhynchus_anatinus	| mOrnAna2	| 1	| 0	| 0	| 0	| X	| X	| X	| X	| dovetail |
| Parambassis_ranga	| fParRan2	| 1	| 0	| 0	| 0	| X	| X	| X	| X	| arima |
| Choloepus_didactylus	| mChoDid1	| 1	| 35	| 35	| 0	| X	| X	| X	| X	|  |
| Nomonyx_dominicus	| bNomDom1	| 1	| 0	| 0	| 0	| O	| X	| O	| X	|  |
| Cottoperca_gobio	| fCotGob2	| 1	| 0	| 0	| 0	| X	| X	| X	| X	| arima |
