#include "DistanceMethods.hpp"
#if defined(_WIN32) || defined(WIN32)
#include "glog/logging.h"
#include <iostream>

std::string RapidNJ(TaxonSet &ts, DistanceMatrix &dm) {

  LOG(ERROR) << "RapidNJ is not supported on Windows" << std::endl;
  return "RapidNJ is not supported on Windows";
}

#else

#include "third_party/rapidNJ/rapidNJ.h"
#include <sstream>

std::string RapidNJ(TaxonSet &ts, DistanceMatrix &dm) {
  bool dmreader_verbose = true;
  int matrixSize = ts.size();
  bool halfMatrix = false;
  std::vector<std::string> sequenceNames;
  distType **matrix;

  for (Taxon t : ts) {
    sequenceNames.push_back(ts[t]);
  }

  matrix = new distType *[matrixSize];
  for (int t = 0; t < matrixSize; t++) {
    matrix[t] = new distType[matrixSize];
    for (int s = 0; s < matrixSize; s++) {
      matrix[t][s] = dm(t, s);
    }
  }

  distMatrixReader *reader = new distMatrixReader(
      dmreader_verbose, matrixSize, halfMatrix, &sequenceNames, matrix);

  ProgressBar pb;
  rapidNJ rnj(reader, matrixSize, false, &pb);
  polytree *tree = rnj.run();

  std::stringstream out;

  tree->serialize_tree(out);

  return out.str();
}

#endif