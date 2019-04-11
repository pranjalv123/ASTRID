# ASTRID-2

This repository is under development, and there may be some changes to the user interface in the future.


ASTRID-2 is a method for estimating species trees from gene trees. 

To build ASTRID-2, clone the git repository, and do the following:

    mkdir build
    cd build
    cmake ../src
    make
   
To run ASTRID-2, you can do

    ASTRID -i <input gene trees> -o <output species tree>
    
to run ASTRID with the BME criterion in FastME, with SPR local search (if the distance matrix is complete) or with BioNJ* (if the distance matrix is incomplete).

Note that if you are trying to run BioNJ*, you must have the PhyDstar.jar file in the same file as the ASTRID executable.

It's possible that FastME with SPR is slow on extremely large datasets, so you can instead do 

    ASTRID -i <input gene trees> -o <output species tree> -n
    
to do NNIs instead of SPRs, or 

    ASTRID -i <input gene trees> -o <output species tree> -f
    
to do no local search. 

On moderately large datasets BIONJ* maybe slow; to run on datasets with high levels of missing data I recommend:

    ASTRID -i <input gene trees> -o <output species tree> -u -n -s
    
which first runs UPGMA on the matrix with missing data, then fills in the missing elements of the distance matrix with distances from the UPGMA tree.
It then repeats the process with the FastME-NNI tree, then gives the final result using FastME-SPR.

This is experimental; I expect the accuracy to be superior to that with BIONJ*, but this is unpublished and undertested.