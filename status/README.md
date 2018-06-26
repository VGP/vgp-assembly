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

### Status as of 2018-06-26

O: Available

X: Not Available

Total num. objects under s3://genomeark/species: 29,103

| species_name	| species_id	| tech_count	| pacbio_subreads	| pacbio_scrubs	| 10x	| bionano_bnx	| bionano_cmap	| hic |
| :---------- | :---------- | :---------- | :---------- | :---------- | :----- | :----- | :----- | :-----  |
| Calypte_anna	| bCalAnn1	| 4	| 189	| 0	| 1	| DLE1,BspQI,BssSI	| BspQI,BssSI	| arima,dovetail,phase |
| Phyllostomus_discolor	| mPhyDis1	| 4	| 43	| 43	| 8	| BspQI,BssSI	| BspQI,BssSI	| arima |
| Rhinatrema_bivittatum	| aRhiBiv1	| 4	| 54	| 54	| 12	| BspQI,BssSI	| BspQI,BssSI	| arima |
| Anabas_testudineus	| fAnaTes1	| 4	| 8	| 8	| 4	| BspQI,BssSI	| BspQI,BssSI	| arima |
| Gopherus_evgoodei	| rGopEvg1	| 4	| 19	| 19	| 4	| DLE1	| DLE1	| arima |
| Rhinolophus_ferrumequinum	| mRhiFer1	| 4	| 25	| 25	| 12	| BspQI,BssSI	| BspQI,BssSI	| phase |
| Strigops_habroptilus	| bStrHab1	| 4	| 126	| 0	| 4	| DLE1	| DLE1	| arima |
| Gouania_willdenowi	| fGouWil2	| 4	| 16	| 16	| 4	| BspQI,BssSI	| BspQI,BssSI	| arima |
| Amblyraja_radiata	| sAmbRad1	| 4	| 50	| 50	| 12	| BspQI,BssSI	| BspQI,BssSI	| arima |
| Lynx_canadensis	| mLynCan4	| 4	| 35	| 35	| 8	| BspQI,BssSI	| BspQI,BssSI	| arima |
| Mastacembelus_armatus	| fMasArm1	| 4	| 8	| 8	| 4	| BspQI,BssSI	| BspQI,BssSI	| arima |
| Taeniopygia_guttata	| bTaeGut2	| 4	| 20	| 20	| 8	| BspQI,BssSI	| X	| arima |
| Taeniopygia_guttata	| bTaeGut1	| 4	| 118	| 0	| 8	| BspQI,BssSI	| X	| arima |
| Geothlypis_trichas	| bGeoTri1	| 4	| 27	| 27	| 12	| BspQI,BssSI	| X	| arima |
| Astatotilapia_calliptera	| fAstCal1	| 3	| 67	| 0	| 4	| BspQI	| BspQI	| X |
| Ornithorhynchus_anatinus	| mOrnAna1	| 3	| 101	| 101	| 4	| DLE1	| DLE1	| X |
| Cottoperca_gobio	| fCotGob3	| 3	| 10	| 10	| 4	| BspQI,BssSI	| BspQI,BssSI	| X |
| Sparus_aurata	| fSpaAur1	| 3	| 10	| 10	| 12	| X	| X	| arima |
| Denticeps_clupeoides	| fDenClu1	| 3	| 8	| 8	| 4	| X	| X	| arima |
| Archocentrus_centrarchus	| fArcCen1	| 3	| 18	| 0	| 4	| X	| X	| phase |
| Takifugu_rubripes	| fTakRub1	| 3	| 4	| 4	| 4	| X	| X	| arima |
| Dermochelys_coriacea	| rDerCor1	| 3	| 17	| 17	| 0	| DLE1	| X	| arima |
| Salmo_trutta	| fSalTru1	| 3	| 26	| 26	| 8	| X	| X	| arima |
| Streptopelia_turtur	| bStrTur1	| 3	| 18	| 18	| 12	| X	| X	| arima |
| Erpetoichthys_calabaricus	| fErpCal1	| 3	| 53	| 53	| 8	| X	| X	| arima |
| Dendrocygna_viduata	| bDenVid1	| 2	| 10	| 10	| 0	| BspQI,BssSI	| X	| X |
| Spatula_cyanoptera	| bSpaCya1	| 2	| 18	| 18	| 0	| BspQI,BssSI	| X	| X |
| Phoenicopterus_ruber	| bPhoRub1	| 2	| 8	| 8	| 0	| DLE1	| X	| X |
| Acanthisitta_chloris	| bAcaChl1	| 2	| 8	| 8	| 0	| DLE1	| X	| X |
| Aquila_chrysaetos	| bAquChr1	| 2	| 16	| 16	| 12	| X	| X	| X |
| Thalassophryne_amazonica	| fThaAma1	| 2	| 20	| 20	| 8	| X	| X	| X |
| Branta_leucopsis	| bBraLeu1	| 2	| 16	| 16	| 0	| BspQI,BssSI	| X	| X |
| Turdus_merula	| bTurMer1	| 2	| 11	| 11	| 0	| BspQI,BssSI,DLE1	| X	| X |
| Zeus_faber	| fZeuFab1	| 2	| 27	| 27	| 4	| X	| X	| X |
| Cariama_cristata	| bCarCri1	| 2	| 8	| 8	| 0	| DLE1	| X	| X |
| Scyliorhinus_canicula	| sScyCan1	| 2	| 34	| 34	| 16	| X	| X	| X |
| Sciurus_vulgaris	| mSciVul1	| 2	| 36	| 36	| 8	| X	| X	| X |
| Nomonyx_dominicus	| bNomDom1	| 2	| 12	| 12	| 0	| BspQI,BssSI	| X	| X |
| Bucorvus_abyssinicus	| bBucAby1	| 1	| 0	| 0	| 0	| DLE1	| X	| X |
| Microcaecilia_unicolor	| aMicUni1	| 1	| 0	| 0	| 12	| X	| X	| X |
| Ornithorhynchus_anatinus	| mOrnAna2	| 1	| 0	| 0	| 0	| X	| X	| dovetail |
| Geotrypetes_seraphini	| aGeoSer1	| 1	| 0	| 0	| 24	| X	| X	| X |
| Parambassis_ranga	| fParRan2	| 1	| 0	| 0	| 0	| X	| X	| arima |
| Choloepus_didactylus	| mChoDid1	| 1	| 35	| 35	| 0	| X	| X	| X |
| Cottoperca_gobio	| fCotGob2	| 1	| 0	| 0	| 0	| X	| X	| arima |
