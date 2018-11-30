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

### Status as of 2018-11-30

X: Not Available

Total num. objects under s3://genomeark/species: 49,759

| genome_name	| species_id	| tech_count	| pacbio_subreads	| pacbio_scrubs	| 10x	| bionano_bnx	| bionano_cmap	| hic 	| assembly |
| :---------- | :---------- | :---------- | :---------- | :---------- | :----- | :----- | :----- | :----- | :----- |
| Calypte_anna	| bCalAnn1	| 5	| 189	| 0	| 1	| BspQI,BssSI,DLE1	| BspQI,BssSI,DLE1	| arima,dovetail,phase	| bCalAnn1.alt.cur.20181019,bCalAnn1.pri.cur.20180926,p1,q2,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,v1.h,v1.p,bCalAnn1.alt.asm.20180817,bCalAnn1.pri.asm.20180817 |
| Rhinatrema_bivittatum	| aRhiBiv1	| 5	| 54	| 54	| 12	| BspQI,BssSI	| BspQI,BssSI	| arima	| aRhiBiv1.alt.cur.20180909,aRhiBiv1.pri.cur.20180909,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,t4.h,t4.p |
| Ornithorhynchus_anatinus	| mOrnAna1	| 5	| 101	| 101	| 4	| DLE1	| DLE1	| phase	| mOrnAna1.alt.cur.20181116,mOrnAna1.pri.cur.20181116.MT,mOrnAna1.pri.cur.20181116,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,t4.h,t4.p,mOrnAna1.alt.asm.20180817,mOrnAna1.pri.asm.20180817 |
| Anabas_testudineus	| fAnaTes1	| 5	| 8	| 8	| 4	| BspQI,BssSI	| BspQI,BssSI	| arima	| v0.h,v0.p,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,v1.h,v1.p,fAnaTes1.alt.asm.20180817,fAnaTes1.pri.asm.20180817 |
| Gopherus_evgoodei	| rGopEvg1	| 5	| 19	| 19	| 4	| DLE1	| DLE1	| arima	| c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,v1.h,v1.p,rGopEvg1.alt.asm.20180817,rGopEvg1.pri.asm.20180817,rGopEvg1.alt.asm.20181126,rGopEvg1.pri.asm.20181126 |
| Amblyraja_radiata	| sAmbRad1	| 5	| 50	| 50	| 12	| BspQI,BssSI	| BspQI,BssSI	| arima	| p1,q2,c1,c2,h1,h2,p1,p2,q2 |
| Lynx_canadensis	| mLynCan4	| 5	| 34	| 34	| 8	| BspQI,BssSI	| BspQI,BssSI	| arima	| c1,c2,s1,s2,s3,s4,t1,t2,t3,t3.p,v0.h,v0.p,q2,s2,mLynCan4.alt.asm.20181119,mLynCan4.pri.asm.20181119 |
| Phyllostomus_discolor	| mPhyDis1	| 5	| 43	| 43	| 8	| BspQI,BssSI	| BspQI,BssSI	| arima	| mPhyDis1.alt.cur.20180907,mPhyDis1.pri.cur.20180907,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,v1,v1.h,v1.p,mPhyDis1.alt.asm.20180817,mPhyDis1.pri.asm.20180817 |
| Taeniopygia_guttata	| bTaeGut2	| 5	| 20	| 20	| 8	| BspQI,BssSI	| BspQI,BssSI	| arima	| bTaeGut2.alt.cur.20181019,bTaeGut2.pri.cur.20181019,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,t4.h,t4.p,v1.MT,bTaeGut2.alt.asm.20180817,bTaeGut2.pri.asm.20180817 |
| Taeniopygia_guttata	| bTaeGut1	| 5	| 118	| 0	| 8	| BspQI,BssSI	| BspQI,BssSI	| arima	| bTaeGut1.alt.cur.20181023,bTaeGut1.pri.cur.20181023,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,v1.h,v1.p,bTaeGut1.alt.asm.20180817,bTaeGut1.pri.asm.20180817 |
| Archocentrus_centrarchus	| fArcCen1	| 5	| 18	| 0	| 4	| DLE1	| DLE1	| phase	| fArcCen1.alt.cur.20180907,fArcCen1.pri.cur.20180907,fArcCen1.pri.cur.20181118,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,t4.h,t4.p,t5.p,fArcCen1.alt.asm.20180817,fArcCen1.pri.asm.20180817,fArcCen1.pri.asm.20181002 |
| Rhinolophus_ferrumequinum	| mRhiFer1	| 5	| 25	| 25	| 12	| BspQI,BssSI	| BspQI,BssSI	| phase	| mRhiFer1.alt.cur.20180907,mRhiFer1.pri.cur.20180907,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,v1,v1.h,v1.p,mRhiFer1.alt.asm.20180817,mRhiFer1.pri.asm.20180817 |
| Strigops_habroptilus	| bStrHab1	| 5	| 126	| 0	| 4	| DLE1	| DLE1	| arima	| bStrHab1.alt.cur.20180907,bStrHab1.alt.cur.20181120.MT,bStrHab1.alt.cur.20181120,bStrHab1.pri.cur.20180907.MT,bStrHab1.pri.cur.20180907,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,v1.MT,v1.h,v1.p,bStrHab1.alt.asm.20180817,bStrHab1.pri.asm.20180817 |
| Gouania_willdenowi	| fGouWil2	| 5	| 16	| 16	| 4	| BspQI,BssSI	| BspQI,BssSI	| arima	| fGouWil2.alt.cur.20181116,fGouWil2.pri.cur.20181116,c1,c2,p1,p2,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,c1,c2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,fGouWil2.alt.asm.20180817,fGouWil2.pri.asm.20180817,fGouWil2.alt.asm.20180830,fGouWil2.pri.asm.20180830 |
| Mastacembelus_armatus	| fMasArm1	| 5	| 8	| 8	| 4	| BspQI,BssSI	| BspQI,BssSI	| arima	| v0.h,v0.p,c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,t4.h,t4.p,fMasArm1.alt.asm.20180817,fMasArm1.pri.asm.20180817 |
| Astatotilapia_calliptera	| fAstCal1	| 4	| 67	| 0	| 4	| BspQI	| BspQI	| X	| v0.h,v0.p,c1,c2,p1,q2,s1,s2,s4,t1,t2,t3,t3.h,t3.p,v1.h,v1.p,fAstCal1.alt.asm.20180817,fAstCal1.pri.asm.20180817 |
| Cottoperca_gobio	| fCotGob3	| 4	| 10	| 10	| 4	| BspQI,BssSI	| BspQI,BssSI	| X	| c1,c2,p1,q2,s1,s2,s3,s4,t1,t2,t3,t3.h,t3.p,fCotGob3.alt.asm.20180817,fCotGob3.pri.asm.20180817 |
| Taeniopygia_guttata	| bTaeGut4	| 2	| 0	| 0	| 0	| X	| X	| illumina	| c |
| Taeniopygia_guttata	| bTaeGut3	| 2	| 0	| 0	| 0	| X	| X	| illumina	| c |
| Denticeps_clupeoides	| fDenClu1	| 4	| 8	| 8	| 4	| DLE1	| DLE1	| arima	| X |
| Geothlypis_trichas	| bGeoTri1	| 4	| 27	| 27	| 12	| BspQI,BssSI	| BspQI,BssSI	| arima	| X |
| Takifugu_rubripes	| fTakRub1	| 4	| 4	| 4	| 4	| DLE1	| DLE1	| arima	| X |
| Salmo_trutta	| fSalTru1	| 4	| 26	| 26	| 8	| DLE1	| DLE1	| arima	| X |
| Parambassis_ranga	| fParRan2	| 4	| 7	| 7	| 4	| DLE1	| DLE1	| arima	| X |
| Erpetoichthys_calabaricus	| fErpCal1	| 4	| 53	| 53	| 8	| DLE1	| DLE1	| arima	| X |
| Phoenicopterus_ruber	| bPhoRub1	| 3	| 8	| 8	| 8	| DLE1	| DLE1	| X	| X |
| Cariama_cristata	| bCarCri1	| 3	| 8	| 8	| 4	| DLE1	| DLE1	| X	| X |
| Scyliorhinus_canicula	| sScyCan1	| 3	| 34	| 34	| 16	| DLE1	| DLE1	| X	| X |
| Sciurus_vulgaris	| mSciVul1	| 3	| 36	| 36	| 8	| DLE1	| DLE1	| X	| X |
| Alca_torda	| bAlcTor1	| 3	| 10	| 10	| 8	| DLE1	| DLE1	| X	| X |
| Bucorvus_abyssinicus	| bBucAby1	| 3	| 13	| 13	| 4	| DLE1	| DLE1	| X	| X |
| Acanthisitta_chloris	| bAcaChl1	| 3	| 8	| 8	| 8	| DLE1	| DLE1	| X	| X |
| Callithrix_jacchus	| mCalJac1	| 3	| 33	| 33	| 4	| DLE1	| DLE1	| X	| X |
| Chiroxiphia_lanceolata	| bChiLan1	| 3	| 10	| 10	| 4	| DLE1	| DLE1	| X	| X |
| Dermochelys_coriacea	| rDerCor1	| 3	| 17	| 17	| 0	| DLE1	| DLE1	| arima	| X |
| Choloepus_didactylus	| mChoDid1	| 3	| 35	| 35	| 4	| DLE1	| DLE1	| X	| X |
| Tachyglossus_aculeatus	| mTacAcu1	| 2	| 0	| 0	| 4	| DLE1	| DLE1	| X	| X |
| Spatula_cyanoptera	| bSpaCya1	| 2	| 18	| 18	| 0	| BspQI,BssSI	| BspQI,BssSI	| X	| X |
| Cygnus_olor	| bCygOlo1	| 2	| 14	| 14	| 0	| DLE1	| DLE1	| X	| X |
| Branta_leucopsis	| bBraLeu1	| 2	| 16	| 16	| 0	| BspQI,BssSI	| BspQI,BssSI	| X	| X |
| Tauraco_erythrolophus	| bTauEry1	| 2	| 9	| 9	| 0	| DLE1	| DLE1	| X	| X |
| Notolabrus_celidotus	| fNotCel1	| 2	| 10	| 10	| 0	| DLE1	| DLE1	| X	| X |
| Nomonyx_dominicus	| bNomDom1	| 2	| 12	| 12	| 0	| BspQI,BssSI	| BspQI,BssSI	| X	| X |
| Dendrocygna_viduata	| bDenVid1	| 2	| 10	| 10	| 0	| BspQI,BssSI	| BspQI,BssSI	| X	| X |
| Sparus_aurata	| fSpaAur1	| 3	| 10	| 10	| 12	| X	| X	| arima	| X |
| Aythya_fuligula	| bAytFul2	| 2	| 7	| 7	| 0	| DLE1	| DLE1	| X	| X |
| Thamnophis_elegans	| rThaEle1	| 2	| 13	| 13	| 0	| DLE1	| DLE1	| X	| X |
| Lacerta_agilis	| rLacAgi1	| 2	| 19	| 19	| 0	| DLE1	| DLE1	| X	| X |
| Antennarius_maculatus	| fAntMac1	| 3	| 24	| 0	| 4	| X	| X	| phase	| X |
| Merops_nubicus	| bMerNub1	| 2	| 10	| 10	| 0	| DLE1	| DLE1	| X	| X |
| Streptopelia_turtur	| bStrTur1	| 3	| 18	| 18	| 12	| X	| X	| arima	| X |
| Rousettus_aegyptiacus	| mRouAeg1	| 1	| 0	| 0	| 0	| DLE1	| DLE1	| X	| X |
| Trichosurus_vulpecula	| mTriVul1	| 2	| 30	| 30	| 4	| X	| X	| X	| X |
| Catharus_ustulatus	| bCatUst1	| 1	| 0	| 0	| 0	| DLE1	| DLE1	| X	| X |
| Sphaeramia_orbicularis	| fSphaOr1	| 2	| 7	| 7	| 12	| X	| X	| X	| X |
| Amphilophus_labiatus	| fAmpLab1	| 1	| 0	| 0	| 0	| DLE1	| DLE1	| X	| X |
| Hippoglossus_hippoglossus	| fHipHip1	| 1	| 0	| 0	| 0	| DLE1	| DLE1	| X	| X |
| Myripristis_murdjan	| fMyrMur1	| 2	| 8	| 8	| 12	| X	| X	| X	| X |
| Zeus_faber	| fZeuFab1	| 2	| 27	| 27	| 4	| X	| X	| X	| X |
| Scleropages_formosus	| fSclFor1	| 2	| 8	| 8	| 12	| X	| X	| X	| X |
| Syngnathus_acus	| fSynAcu1	| 2	| 3	| 3	| 4	| X	| X	| X	| X |
| Aquila_chrysaetos	| bAquChr1	| 2	| 16	| 16	| 12	| X	| X	| X	| X |
| Salarias_fasciatus	| fSalaFa1	| 2	| 13	| 13	| 4	| X	| X	| X	| X |
| Thalassophryne_amazonica	| fThaAma1	| 2	| 20	| 20	| 8	| X	| X	| X	| X |
| Turdus_merula	| bTurMer1	| 2	| 11	| 11	| 0	| BspQI,BssSI,DLE1	| X	| X	| X |
| Microcaecilia_unicolor	| aMicUni1	| 2	| 36	| 36	| 12	| X	| X	| X	| X |
| Echeneis_naucrates	| fEcheNa1	| 2	| 5	| 5	| 4	| X	| X	| X	| X |
| Lutra_lutra	| mLutLut1	| 2	| 16	| 16	| 8	| X	| X	| X	| X |
| Geotrypetes_seraphini	| aGeoSer1	| 2	| 6	| 6	| 24	| X	| X	| X	| X |
| Dendropsophus_ebraccatus	| aDenEbr1	| 1	| 0	| 0	| 0	| DLE1	| DLE1	| X	| X |
| Sylvia_atricapilla	| bSylAtr1	| 2	| 8	| 8	| 8	| X	| X	| X	| X |
| Balaenoptera_musculus	| mBalMus1	| 1	| 52	| 52	| 0	| X	| X	| X	| X |
| Ornithorhynchus_anatinus	| mOrnAna2	| 1	| 0	| 0	| 0	| X	| X	| dovetail	| X |
| Cottoperca_gobio	| fCotGob2	| 1	| 0	| 0	| 0	| X	| X	| arima	| X |
| Xenentodon_cancila	| fXenCan1	| 1	| 14	| 0	| 0	| X	| X	| X	| X |
| Anableps_anableps	| fAnaAna1	| 1	| 13	| 0	| 0	| X	| X	| X	| X |

