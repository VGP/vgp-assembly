## Overview

This script is to trim off large stretches of N bases and to remove fake cut sites that sometimes get inserted at the end of bionano hybrid scaffolds.

Tested on: Python3.4.0 (may work with any Python3), requires BioPython

Step 1. First you run `remove_fake_cut_sites_DNAnexus.py` to generate a FASTA 
file from which BioNano cut sites in runs of Ns have been removed.

Step 2. Then you run that FASTA file through `trim_Ns_DNAnexus.py` to obtain a 
list of regions with Ns. (`trim_Ns_DNAnexus.pl` is an alternative for this stage)

Step 3. Finally, you combine the FASTA file from step 1 with the list from step 2 to obtain a trimmed FASTA file.


## Details of each script


### remove_fake_cut_sites_DNAnexus.py

This script removes the spurious BioNano cut sites that are inserted 
into gaps in some assemblies, replacing them with Ns.

```
Usage: python3 remove_fake_cut_sites_DNAnexus.py <input.fa> <output.fa> <output.log>
```

Where:

* `input.fa` is the input FASTA file

* `output.fa` is the output with any cut sites in the middle of a run of 
Ns replaced with Ns.

* `output.log` lists the number of cut sites replaced in each scaffold.

Note that this only replaces cut sites belonging to a short list given 
in the script (and their reverse complements); other cut sites will not 
be replaced. However, the log file records the presence of any runs of 
sequence between 1 and 10 bp that are flanked by Ns which *don't* get 
replaced, which should alert one to the presence of fake cut sites.



### trim_Ns_DNAnexus.py

This identifies runs of Ns near the ends of sequences (despite the name 
of the script, it only identifies these regions; it does not actually 
trim them). More specifically, it identifies runs of Ns that are 
*directly* at the ends of scaffolds, but it also removes runs of Ns that 
are very near the ends of scaffolds (it looks at the 5 kb window at each 
end of the scaffold, and if more than 3 kb of that is Ns, it slides the 
window towards the center of the sequence until it finds a 5kb window 
where less than 3kb is Ns, and then trims everything up to that point 
other than the final run of non-Ns). Furthermore, it flags any scaffold 
for removal if it has fewer than 100 non-N bases. The output format is a 
bit idiosyncratic and includes some redundancies where sequences fit 
several of the categories given above, but the clipping script below 
takes this into account.

```
Usage: python3 trim_Ns_DNAnexus.py <input.fa> <output.list>
```

Where:

* `input.fa` is the input FASTA file

* `output.list` is a record of the locations for trimming

### trim_Ns_DNAnexus.pl

This is a Perl alternative for the Python script given above. It's 
actually our earlier version of the script. The output should be 
identical. The Perl script is generally slower, and under some 
circumstances can be *much* slower, although it's generally 
satisfactory. The reason I'm including it as an option is that the 
Python script is a bit baroque (in order to be faster than the Perl 
script but obey the same logic), whereas the Perl script is quite 
straightforward. I discovered a bug in the Python script yesterday that 
triggered under rare circumstances; that's fixed now, but it reminded me 
that the non-intuitive nature of the Python algorithm meant it might 
conceivably harbour other obscure bugs, and made me feel I should offer 
the Perl script as an alternative if you want to play things safer.

The Perl script takes the same arguments and gives the same output as 
the Python script.

### clip_regions_DNAnexus.py

This exists as a separate script from those above because in our 
pipeline it represents the final stage of contamination-checking 
including not only the trim_Ns script, but also other forms of 
contamination-checking. For that reason, a good deal of this script is 
concerned with contingencies other than trimming Ns which won't be 
relevant here.

Usage:

```
python3 clip_regions_DNAnexus.py <input.fa> <input.list> <output.fa>
```

Where:

* `input.fa` is the file to be trimmed (in this case, the output from 
remove_fake_cut_sites_DNAnexus.py)

* `input.list` is the list of regions to be trimmed (the output from 
trim_Ns_DNAnexus.py)

* `output.fa` is the name for the trimmed FASTA file to be created.

*This README file was generated based on the email conversation with James Torrance. March 13, 2019*

