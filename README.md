# Recombination Simulation

This repository contains scripts to create a multiple-sequence alignment of recombinant genomes. It was designed to test the recombination module in [UShER](https://usher-wiki.readthedocs.io/en/latest/), but may be used for other purposes as well.

#### Summary:

- *makeInternalNodeMSA.py* generates a MSA of the internal nodes from the initial tree with at least 10 descendants.
--- test
- *makeRandomRecombinants.py* generates a MSA of recombinant samples from the internal node sequences.
- *parseRecombinationResults.py* reads in results from findRecombination, as well as the log and msa corresponding to the recombinant sequences, and prints the recombinant samples whose breakpoint was predicted incorrectly.

#### Example pipeline:

```
matUtils extract -i input.pb -A samples.tsv -L num_leaves.tsv
python makeInternalNodeMSA.py -s samples.tsv -l num_leaves.tsv -r wuhan.ref.fa
python makeRandomRecombinants.py -b 1 -s 100 -c 10 -d 10 -f internal_nodes.msa.fa -r wuhan.ref.fa
faToVcf recombination_1.msa.fa recombination_1.msa.vcf
usher -i input.pb -v recombination_1.msa.vcf -o recombination_1.msa.pb
findRecombination -i recombination_1.msa.pb
python parseRecombinationResults.py -f recombination_1.msa.fa -l recombination_1.log -d descendants.tsv -r recombination.tsv
```
