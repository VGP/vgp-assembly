# mitoVGP 2.2
This repository contains scripts used to generate mitochondrial sequences for the <a href="http://www.vertebrategenomesproject.org">Vertebrate Genomes Project</a>.

<b>Software and Data Use Policy</b>

mitoVGP is distributed under the <a href="LICENSE.txt">BSD 3-Clause License</a>.

VGP samples and data come from a variety of sources. To support fair and productive use of this data, please abide by the <a href="https://genome10k.soe.ucsc.edu/data-use-policies/">Data Use Policy</a> and contact us with any questions.

<b>Content Description:</b>

- canu-1.8.Linux-amd64.tar.xz - the popular long read assembler employed in the pipeline

- mitoVGP_conda_env_pacbio.yml - conda environment containing all software required to run the pipeline with Pacbio data on Linux

- mitoVGP_conda_env_ONT.yml - conda environment containing all software required to run the pipeline with ONT data on Linux

- mitoVGP - the pipeline

- scripts/ - the intermediate scripts required by mitoVGP

<b>Quick Start</b>

mitoVGP is available for Linux64 and requires <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation">Conda</a>. To install and run follow these instructions:

```
git clone https://github.com/gf777/mitoVGP.git #clone this git repository
cd mitoVGP #get into mitoVGP folder

tar -xvf canu-1.8.Linux-amd64.tar.xz #install canu assembler
rm canu-1.8.Linux-amd64.tar.xz

#create the mitoVGP pipeline software environment
#please note: Pacbio software only runs on Python 2, while ONT software requires Python 3,
#therefore two different environments must be set depending on data type.
#Pacbio:
conda env create -f mitoVGP_conda_env_pacbio.yml
#ONT:
conda env create -f mitoVGP_conda_env_ONT.yml

conda activate mitoVGP_pacbio #activate mitoVGP conda environment, use mitoVGP_ONT

#run mitoVGP pipeline using 24 cores (example with M. armatus, Pacbio data)
./mitoVGP -a pacbio -s Mastacembelus_armatus -i fMasArm1 -r mtDNA_Mastacembelus_armatus.fasta -t 24 -b variantCaller
```

For additional options and specifications you can type:
```
./mitoVGP -h
```

Please note that depending on your Pacbio chemistry you will need to define a different polishing tool.
For chemistry 2.0 (default):
```
./mitoVGP -b gcpp
```
For chemistry lower than 2.0 add:
```
./mitoVGP -b variantCaller
```
For very old RSII chemistries you may also want to align reads using blasr:
```
./mitoVGP -b variantCaller -m blasr
```


<b> Pipeline workflow </b>

An existing reference from closely to distantly related species is used to identify mito-like reads in pacbio/ONT WGS data, which are then employed in <i>de novo</i> genome assembly. The assembly is further polished using both long and short read data, and linearized to start with the conventional Phenylalanine tRNA sequence.

<p align="center">
	<img src="mitoVGP_pipeline_Rockefeller_v.2.2.png" />
</p>

VGP mitogenomes assembled using mitoVGP pipeline can be found on <a href="https://vgp.github.io/genomeark/">GenomeArk</a> and include:

<b>Pacbio</b><br/>
<i>
Acanthisitta chloris<br/>
Acipenser ruthenus<br/>
Alca torda<br/>
Amblyraja radiata<br/>
Anabas testudineus<br/>
Anableps anableps<br/>
Antennarius maculatus<br/>
Aquila chrysaetos<br/>
Archocentrus centrarchus<br/>
Arvicanthis niloticus<br/>
Astatotilapia calliptera<br/>
Asterias rubens<br/>
Balaenoptera musculus<br/>
Bos taurus<br/>
Bufo bufo<br/>
Callithrix jacchus<br/>
Calypte anna<br/>
Carcharodon carcharias<br/>
Cariama cristata<br/>
Catharus ustulatus<br/>
Chelmon rostratus<br/>
Choloepus didactylus<br/>
Ciconia maguari<br/>
Corvus moneduloides<br/>
Cyclopterus lumpus<br/>
Cygnus olor<br/>
Danio rerio<br/>
Dendropsophus ebraccatus<br/>
Denticeps clupeoides<br/>
Dermochelys coriacea<br/>
Dryobates pubescens<br/>
Echeneis naucrates<br/>
Electrophorus electricus<br/>
Erpetoichthys calabaricus<br/>
Esox lucius<br/>
Falco naumanni<br/>
Falco rusticolus<br/>
Gallus gallus<br/>
Geotrypetes seraphini<br/>
Gopherus evgoodei<br/>
Gouania willdenowi<br/>
Hemiprocne comata<br/>
Hippoglossus hippoglossus<br/>
Homo sapiens<br/>
Lacerta agilis<br/>
Lemur catta<br/>
Lutra lutra<br/>
Lynx canadensis<br/>
Mastacembelus armatus<br/>
Megalops cyprinoides<br/>
Melanotaenia boesemani<br/>
Melopsittacus undulatus<br/>
Merops nubicus<br/>
Microcaecilia unicolor<br/>
Mustela erminea<br/>
Myotis myotis<br/>
Notolabrus celidotus<br/>
Nyctibius grandis<br/>
Ornithorhynchus anatinus<br/>
Pan troglodytes<br/>
Periophthalmus magnuspinnatus<br/>
Phocoena sinus<br/>
Phyllostomus discolor<br/>
Pipistrellus kuhlii<br/>
Pipistrellus pipistrellus<br/>
Pluvialis apricaria<br/>
Pristis pectinata<br/>
Pterocles gutturalis<br/>
Pygocentrus nattereri<br/>
Rana temporaria<br/>
Rattus norvegicus<br/>
Rhinatrema bivittatum<br/>
Rhinolophus ferrumequinum<br/>
Salmo trutta<br/>
Scatophagus argus <br/>
Sciurus vulgaris<br/>
Sebastes umbrosus<br/>
Sparus aurata<br/>
Sterna hirundo<br/>
Strigops habroptilus<br/>
Syngnathus acus<br/>
Tachyglossus aculeatus<br/>
Taeniopygia guttata<br/>
Taeniopygia guttata<br/>
Takifugu rubripes<br/>
Trachurus trachurus<br/>
Trichosurus vulpecula<br/>
Tursiops truncatus<br/>
Xenentodon cancila<br/>
Zalophus californianus<br/>
Zeus faber<br/>
</i>

<br/>

<b>Nanopore</b><br/>
<i>
Notolabrus celidotus<br/>
</i>