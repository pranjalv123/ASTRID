#include "DistanceMethods.hpp"

extern "C" {
#include "fastme/fastme.h"
}

#include <sstream>
#include <newick.hpp>

typedef set mySet;

string FastME (TaxonSet& ts, DistanceMatrix& dm, int nni, int spr)
{
	int size = ts.size();
	int numSpecies = size;
	double **D, **A;

	A = initDoubleMatrix (2*numSpecies-2);
	D = initDoubleMatrix (2*numSpecies-2);
	fillZeroMatrix (&A, 2*numSpecies-2);


	for (int i=0; i<size; i++)
	{
		D[i] = (double *) mCalloc (size, sizeof (double));
	}

	for (int i = 0; i < size; i++) {
	  for (int j = 0; j < size; j++) {
	    if (!dm.isMasked(i,j))
	      D[i][j] = dm(i,j);
	  }
	}


	Options options;
	Set_Defaults_Input (&options);
	options.method = TaxAddBAL;
	options.use_SPR = spr;
	options.use_NNI = nni;
	options.NNI    = BALNNI;

	mySet species;
	species.firstNode = 0;
	species.secondNode = 0;

	for (Taxon t : ts) {

	  stringstream ss;
	  ss << t;
	  node*  v = makeNode (ss.str().c_str(), -1);
	  v->index2 = t;
	  addToSet(v, &species);
	}


	tree* t = ComputeTree (&options, D, A, &species, numSpecies, 8);
	int nnicount;
	int sprcount;
	t = ImproveTree(&options, t, D, A, &nnicount, &sprcount, options.fpO_stat_file);

	char* tree_output = new char[size << 10];
	tree_output[0] = '\0';

	NewickPrintTreeStr (t, tree_output, 2);
	return unmap_newick_names(string(tree_output), ts);
}
