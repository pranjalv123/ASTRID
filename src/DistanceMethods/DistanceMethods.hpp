#ifndef ASTRID_DISTANCE_METHODS__
#define ASTRID_DISTANCE_METHODS__


#include "phylokit/DistanceMatrix.hpp"
#include <queue>
#include <sstream>


std::string FastME (TaxonSet& ts, DistanceMatrix& dm, int nni, int spr);
std::string UPGMA(TaxonSet& ts, DistanceMatrix& dm);

#endif
