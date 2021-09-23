# ASTRID-2

This repository is under development, and there may be some changes to the user interface in the future.


ASTRID-2 is a method for estimating species trees from gene trees. 

To build ASTRID-2, clone the git repository, and do the following:

    bazel build //src:ASTRID -c opt

You can install the bazel build system [here](http://bazel.build). Once the build is complete, your executable will be located at bazel-bin/src/ASTRID.
   
To run ASTRID-2, you can do

    ASTRID -i <input gene trees> -o <output species tree>
    
to run ASTRID with the BME criterion in FastME, with SPR local search (if the distance matrix is complete), first completing the tree with UPGMA* if needed.

More command-line options:

    -i --input <filename>       Specify input file. Should contain one tree per line, in Newick format.
    -o --output <filename>      Specify output file (default is infile.astrid). Will contain one output tree, in Newick format. The tree may have branch lengths, depending on the distance method used, but these are not meaningful
    --cache <filename>          Output the distance matrix. This is useful if you want to run a different distance method.  

    [Distance method selection]
    -u       Use UPGMA as distance method for tree estimation
    -f       Use FastME with no local search as distance method for tree estimation
    -n       Use FastME with NNIs for local search as distance method for tree estimation
    -s       Use FastME with NNIs and SPRs for local search as distance method for tree estimation
    --bionj          Use BioNJ* with as distance method for tree estimation (make sure PhyDstar.jar is in the same folder as the ASTRID executable)
    --auto   [default] Automatically choose between --bionj and -s depending on if the distance matrix is missing taxa

    [Multiple individuals]
    -a --multind     [experimental] Specify mapping file for multiple-individual datasets. Can be in one of three formats, which will be auto-detected:

            Format 1:
                    species1:indiv1,indiv2,indiv3
                    species2:,indiv4,indiv5
                    ...

            Format 2:
                    species1 indiv1 indiv2 indiv3
                    species2 indiv4 indiv5
                    ...
            Format 3:
                    indiv1 species1
                    indiv2 species1
                    indiv3 species1
                    indiv4 species2
                    indiv5 species2
                    ...



Note that if you are trying to run BioNJ*, you need to build a different target. 
THIS IS NOT RECOMMENDED AND IS ONLY INCLUDED FOR HISTORICAL REASONS.

    bazel build //src:ASTRID-phydstar -c opt

It's possible that FastME with SPR is slow on extremely large datasets, so you can instead do 

    ASTRID -i <input gene trees> -o <output species tree> -n
    
to do NNIs instead of SPRs, or 

    ASTRID -i <input gene trees> -o <output species tree> -f
    
to do no local search. 

On moderately large datasets BIONJ* maybe slow; to run on datasets with high levels of missing data I recommend:

    ASTRID -i <input gene trees> -o <output species tree> -u -n -s
    
which first runs UPGMA on the matrix with missing data, then fills in the missing elements of the distance matrix with distances from the UPGMA tree.
It then repeats the process with the FastME-NNI tree, then gives the final result using FastME-SPR.

ASTRID now supports multi-individual gene trees. The matrices corresponding to each gene tree are "collapsed" by averaging the distances between individuals of different species. This technique is described in more detail in ["Species Tree Inference from Gene Splits by Unrooted STAR Methods" by Allman, Degnan, and Rhodes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388605/). 
