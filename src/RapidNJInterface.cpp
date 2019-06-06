#include "DistanceMethods.hpp"
#include "rapidNJ/src/rapidNJ.h"
#include <sstream>

string RapidNJ (TaxonSet& ts, DistanceMatrix& dm)
{
  bool dmreader_verbose = true;
  int matrixSize = ts.size();
  bool halfMatrix = false;
  vector<string> sequenceNames;
  distType** matrix;

  for (Taxon t : ts) {
    sequenceNames.push_back(ts[t]);
  }

  matrix = new distType*[matrixSize];
  for (int t = 0; t < matrixSize; t++) {
    matrix[t] = new distType[matrixSize];
    for (int s = 0; s < matrixSize; s++) {      
      matrix[t][s] = dm(t,s);
    }
  }



  
  
  distMatrixReader* reader = new distMatrixReader(dmreader_verbose, matrixSize, halfMatrix, &sequenceNames, matrix);


  ProgressBar pb;
  rapidNJ rnj(reader, matrixSize, false, &pb);
  polytree* tree = rnj.run();
  

  std::stringstream out;

  tree->serialize_tree(out);

  return out.str();
  
}

