# VGP specimen naming scheme

Individual specimens to be seqeunced as part of the Vertebrate Genomes Project (VGP) will be assigned a VGP ID according to the following scheme.
The ID will take the form:
```
[abfmrs]AbcXyz{#}
```
where

  * The one letter prefix `[abfmrs]` corresponds to one of:

    | prefix | class                | 
    |:-------|:---------------------|
    | a      | amphibians           | 
    | b      | birds                | 
    | f      | fishes               | 
    | m      | mammal               | 
    | r      | reptiles             | 
    | s      | sharks and relatives | 

  * The six letter combination `AbcXyz` is a species/strain designator.
    In most cases, this will be `GenSpe` for **Gen**us/**Spe**cies but that is not required (see below for how to resolve clashes).

  * `{#}` is an incremental number per individual specimen from the same species.

For each species in the VGP ordinal project, the 7-letter prefix `[abfmrs]AbcXyz` will be pre-assigned to avoid conflicts.
Across all species, there will be clashes for this 7-letter prefix.
There are a few options for dealing with these:

1. Allow the clashes and with individuals/species disambiguated by the final incremental number.
2. Allow variation within the six letter species designator, e.g. a 2-4 split (`GeSpec`), or modified capitalisation (`GENSpe`)

The EBI are in the process of setting up a registry where VGP IDs can be assigned and avoid individual IDs clashing between centres.


## Examples

| VGP ID   | Species and common name                      |
|:---------|:---------------------------------------------|
| fGouWil2 | Gouania willdenowi; blunt-snouted clingfish  | 
| mLemCat1 | Lemur catta; ring-tailed lemur               | 
| aRhiDar3 | Rhinoderma darwinii; Darwin's frog           | 
| bCalAnn1 | Calypte anna; Anna's hummingbird             | 
| rDerCor1 | Dermochelys coriacea; leatherback sea turtle | 
| sCarCar1 | Carcharodon carcharias; great-white shark    | 


## Tissue samples

For a single individual, there may be multiple tissue samples used for transcriptome sequencing.
The proposed scheme to distinguish these samples is:
```
[abfmrs]AbcXyz{#}.tissue{#}
```
where `tissue` should come from an agreed list of terms (to be decided). Examples: `fGouWil2.brain1`, `fGouWil2.eye2`.

If the tissue used for transcriptome sequencing is from a different indiviual than the one sequenced to produce the assembly, then an new individual VGP ID should registered.


## Biosamples

Having assigned VGP IDs, a [BioSamples](https://www.ebi.ac.uk/biosamples/) accession ID should also be generated for the individual.
Agreed metadata (to be decided) should be attached to the BioSamples entries.
Tissue samples should be assigned metadata based on an agreed ontology such as [Uberon](https://www.ebi.ac.uk/ols/ontologies/uberon) and should used the `Derived from` linking facility in BioSamples to indicate the individual source of that tissue sample.


## Extension beyond vertebrates

If this scheme were to extend beyond vertebrates in the VGP, the below is a proposal which would use all the letters of the alphabet to cover the Tree of Life.
This is meant as a pragmatic division rather then a strict taxonomic one.


| prefix | class                        | count  | group            | notes                                                                                     | 
|:-------|:-----------------------------|:-------|:-----------------|:------------------------------------------------------------------------------------------| 
| a      | amphibians                   | 6439   | chordates        |                                                                                           | 
| b      | birds                        | 10301  | chordates        |                                                                                           | 
| c      | non-vascular plants          | 14222  | plants           |                                                                                           | 
| d      | dicotyledons                 | 200000 | plants           | not monophyletic                                                                          | 
| e      | echinoderm                   | 6753   | other animals    |                                                                                           | 
| f      | fishes                       | 31862  | chordates        | lobe-finned and ray finned = Osteichthyes = Teleostomi (excluding tetrapods)              | 
| g      | fungi                        | 123126 | other eukaryotes |                                                                                           | 
| h      | platyhelminths               | 9164   | other animals    |                                                                                           | 
| i      | insects                      | 795000 | other animals    |                                                                                           | 
| j      | jellyfish and other cnidaria | 9747   | other animals    |                                                                                           | 
| k      | other chordates              | 1926   | chordates        | cephalochordates, urochordates (tunicates), jawless fish; not monophyletic                | 
| l      | monocotyledons (lilies etc.) | 51595  | plants           | 'l' for lily                                                                              | 
| m      | mammals                      | 4863   | chordates        |                                                                                           | 
| n      | nematodes                    | 3455   | other animals    |                                                                                           | 
| o      | sponges                      | 8499   | other animals    |                                                                                           | 
| p      | protists                     | 12695  | other eukaryotes | defined here as eukaryotes not animals or plants or fungi; not monophyletic               | 
| q      | other arthropods             | 120000 | other animals    | not insects; not monophyletic                                                             | 
| r      | reptiles                     | 9789   | chordates        | excluding birds                                                                           | 
| s      | sharks and relatives         | 1149   | chordates        | Chondricthyes = Elasmobranchs and Chimaeras                                               | 
| t      | other animal phyla           | 165    | other animals    |                                                                                           | 
| u      | algae                        | 2056   | plants           | not monophyletic                                                                          | 
| v      | other vascular plants        | 66717  | plants           | ferns, cycads, conifers, gingko etc.; not monophyletic                                    | 
| w      | annelids (worms)             | 12738  | other animals    |                                                                                           | 
| x      | molluscs                     | 41646  | other animals    | the "scs" in "moluscs" sounds a bit like it contains an 'x'                               | 
| y      | bacteria                     | 6468   | prokaryotes      |                                                                                           | 
| z      | archea                       | 281    | prokaryotes      | mosses, liverworts, hornworts; not monophyletic                                           | 
| -      | viruses                      |        |                  | '-' for missing                                                                           | 


Equivalently, presented by group:

<table>
  <tr>
    <th>group</th>
    <th>prefix</th>
    <th>class</th>
  </tr>
  <tr>
    <td rowspan="7">chordates (including vertebrates)</td>
    <td>m</td>
    <td>mammals</td>
  </tr>
  <tr>
    <td>b</td>
    <td>birds</td>
  </tr>
  <tr>
    <td>r</td>
    <td>reptiles</td>
  </tr>
  <tr>
    <td>a</td>
    <td>amphibians</td>
  </tr>
  <tr>
    <td>f</td>
    <td>fishes</td>
  </tr>
  <tr>
    <td>s</td>
    <td>sharks</td>
  </tr>
  <tr>
    <td>k</td>
    <td>other chordates</td>
  </tr>
  <tr>
    <td rowspan="10">other animals</td>
    <td>e</td>
    <td>echinoderms</td>
  </tr>
  <tr>
    <td>x</td>
    <td>molluscs</td>
  </tr>
  <tr>
    <td>i</td>
    <td>insects</td>
  </tr>
  <tr>
    <td>q</td>
    <td>other arthropods</td>
  </tr>
  <tr>
    <td>w</td>
    <td>annelids (worms)</td>
  </tr>
  <tr>
    <td>n</td>
    <td>nematodes</td>
  </tr>
  <tr>
    <td>h</td>
    <td>platyhelminths</td>
  </tr>
  <tr>
    <td>j</td>
    <td>jellyfish and other cnidaria</td>
  </tr>
  <tr>
    <td>o</td>
    <td>sponges</td>
  </tr>
  <tr>
    <td>t</td>
    <td>other animal phyla</td>
  </tr>
  <tr>
    <td rowspan="5">plants</td>
    <td>d</td>
    <td>dicotyledons</td>
  </tr>
  <tr>
    <td>l</td>
    <td>monocotyledons</td>
  </tr>
  <tr>
    <td>v</td>
    <td>other vascular plants</td>
  </tr>
  <tr>
    <td>c</td>
    <td>non-vascular plants</td>
  </tr>
  <tr>
    <td>u</td>
    <td>algae</td>
  </tr>
  <tr>
    <td rowspan="2">other eukaryotes</td>
    <td>g</td>
    <td>fungi</td>
  </tr>
  <tr>
    <td>p</td>
    <td>protists</td>
  </tr>
  <tr>
    <td rowspan="2">prokaryotes</td>
    <td>y</td>
    <td>bacteria</td>
  </tr>
  <tr>
    <td>z</td>
    <td>archaea</td>
  </tr>
</table>
