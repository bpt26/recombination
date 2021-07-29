# Recombination Filtering

This directory contains scripts used for filtering false positives by various methods. Note that some scripts contain hard-coded directories that would require editing to run in the user's environment.

We first ran findRecombination on each node in the protobuf, using a heavily parallelized set of commands following the formula: `findRecombination -i optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.pb -n 2 -S <start_node> -E <end_node> -d <node_dir> -T 4 2> <node_log>`, then combined our results using the following commands:

`cat recombination_*.tsv > catRecombination.txt`  
`cat descendants_*.tsv > catDescendants.txt`  

We then ran `python combineAndGetPVals.py` to get only the largest parsimony score improvement for each recombinant node, as well as combine adjacent and overlapping predicted breakpoint intervals, and create a new table including phylogenetic p-values and a list of descendants for each node. Following this, we used `python getAllNodes.py` to retrieve a list of all relevant nodes and their parents, `matUtils extract -i ../optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.pb -s allRelevantNodes.txt -v allRelevantNodes.vcf -T 10` to extract a .vcf containing all of these nodes, and `python getABABA.py` to extract the informative sites from the .vcf for each putative trio.

From there, we used `checkmutant.py` in conjunction with text files output from previous steps to check each trio for sequence quality. We used `awk '$19 == "False"' report_7_9.txt | awk '$14 == "False"' | awk '$11 == "False"' > final_report.txt` to produce a final superset of trios, and `python finishFiltration.py` to filter by p-values, remove recombinants with highly clustered mutations, redundant trios, and circular trios (in which two or more recombinant nodes are parents of each other). The output of this final script produces the set of trios presented in our manuscript.


