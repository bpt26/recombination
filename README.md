# Recombination Simulation

This repository contains scripts to create a multiple-sequence alignment of recombinant genomes. It was designed to test the recombination module in [UShER](https://usher-wiki.readthedocs.io/en/latest/), but may be used for other purposes as well.

## Summary:

### makeInternalNodesMSA.py
generates a MSA of the internal nodes from the initial tree with at least 10 descendants.  
Options:
- -s (--samples): Sample mutation paths file, generated by e.g. *matUtils extract -i <input.pb> -A <sample-paths.tsv>*.  
- -l (--leaves): Leaves per node, generated by e.g. *matUtils extract -i <input.pb -L <num-leaves.tsv>*.  
- -r (--reference): Reference genome, used to generate VCF. Default = 'wuhan.ref.fa' (in this repository, corresponds to [the reference genome of SARS-CoV-2](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) with [problematic sites](https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf) masked.).

### makeRandomRecombinants.py
generates a MSA of recombinant samples from the internal node sequences.  
Options:
- -b (--breakpoints): Number of breakpoints that each recombinant sample will have. Must be 1 or 2 (Default = 1).  
- -s (--samples): Number of recombinant samples to create (Default = 100).  
- -c (--copies): Number of identical copies to make for each recombinant sample (Default = 10).  
- -d (--differences): Minimum mutational distance for acceptor/donor samples (Default = 10).  
- -m (--commonMutations): Number of mutations to add to each copy, shared by all in a set. (Default = 0).  
- -M (--randomMutations): Number of mutations to add to each copy, randomly chosen for each copy. (Default = 0).  
- -f (--fasta): Fasta file containing sequences for acceptor/donor samples. [REQUIRED]  
- -r (--reference): Fasta file containing reference genome for use in creating VCF. (Default = 'wuhan.ref.fa').  


### parseRecombinationResults.py
reads in results from findRecombination, as well as the log and msa corresponding to the recombinant sequences, and prints the recombinant samples whose breakpoint was predicted incorrectly.  
Options:
- -f (--fasta): Fasta file containing sequences for recombinant genomes. [REQUIRED]  
- -l (--log): Log file containing recombinant genomes. Format: recombinantSample sample1 sample2 bp1 (bp2).... [REQUIRED]  
- -d (--descendants): Descendants file output by findRecombination. (Default = 'descendants.tsv').  
- -r (--recombination): Recombination file output by findRecombination. (Default = 'recombination.tsv').  
- -b (--breakpoints): Number of breakpoints that each recombinant sample has. Must be 1 or 2 (Default = 1).  

## Example pipeline:

```
matUtils extract -i input.pb -S samples.tsv -L num_leaves.tsv
python makeInternalNodesMSA.py -s samples.tsv -l num_leaves.tsv -r wuhan.ref.fa
python makeRandomRecombinants.py -b 1 -s 100 -c 10 -d 10 -f internal_nodes.msa.fa -r wuhan.ref.fa
faToVcf recombination_1.msa.fa recombination_1.msa.vcf
usher -i input.pb -v recombination_1.msa.vcf -o recombination_1.msa.pb
findRecombination -i recombination_1.msa.pb
python parseRecombinationResults.py -f recombination_1.msa.fa -l recombination_1.log -d descendants.tsv -r recombination.tsv -b 1
```
