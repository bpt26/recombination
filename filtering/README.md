# Recombination Filtering

This directory contains scripts used for filtering false positives by various methods. Note that some scripts contain hard-coded directories that would require editing to run in the user's environment.

We first ran RIPPLES on sets of 10 nodes in the protobuf, in parallel, in a set of commands following the formula: 


```
mkdir <node>  
ripples -i optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.pb -n 2 -S <start_node> -E <end_node> -d <node_dir> -T 4 2> <node_log>  
mv node_dir/recombination.tsv FINAL_RESULTS/recombination_<node>.tsv
mv node_dir/descendants.tsv FINAL_RESULTS/descendants_<node>.tsv
```

This results in separate, non-overlapping files containing the recombinants and descendants for each set of nodes, without worry of different commands writing to the same file. To combine all of the `recombination.tsv` and `descendants.tsv` files into a single usable output after RIPPLES has searched the entire tree:

```
cd FINAL_RESULTS/  
cat recombination_*.tsv > catRecombination.txt  
cat descendants_*.tsv > catDescendants.txt   
```

From here, we begin post-processing:

```
# Merge adjacent/overlapping breakpoint predictions, remove superfluous headers, add p-values from null distributions:  
python combineAndGetPVals.py 

# Get node information and recombination-informative sites for every putative recombinant, donor, and acceptor:  
matUtils extract -i optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.pb -S sample_paths.txt  
# ^Note: use UShER commit a6f65ade7a6606ef75902ee290585c6db23aeba6 for the above step due to recent formatting changes to -S output file
python getAllNodes.py  
awk '{print "node_"$0}' allRelevantNodes.txt > allRelevantNodeNames.txt  
matUtils extract -i optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.pb -s allRelevantNodes.txt -v allRelevantNodes.vcf -T 10  
python getABABA.py  
python getDescendants.py   
python getSampleInfo.py  

# Check each trio for sequence quality and filter dubious trios (Note: This requires access to GISAID and curation/alignment of raw reads)   
python checkmutant.py <trio-number> -r  
awk '$19 == "False"' report_7_9.txt | awk '$14 == "False"' | awk '$11 == "False"' > final_report.txt  

# Get 3seq p-values and further filter results according to sequence quality, p-values, clustering of mutations, and redundant trios:  
python makeMNK.py  
python checkClusters.py  
awk '$21 <= .20 {print}' combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters.txt > combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters3SeqP02.txt  
python doNewTieBreakers.py  
python removeRedundant.py  
```

The final output of this process should yield all information about the final set of recombinant trios, heavily filtered. We are currently working on making this process more generalizable to other protobufs, as well as condensing into fewer programs.    


