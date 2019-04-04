#ifndef ASTRID_DISTANCE_METHODS__
#define ASTRID_DISTANCE_METHODS__

extern "C" {
#include "fastme/fastme.h"
}
typedef set mySet;
#include <DistanceMatrix.hpp>
#include <queue>
#include <sstream>


string FastME (TaxonSet& ts, DistanceMatrix& dm, int nni, int spr);


string BioNJStar(TaxonSet& ts, DistanceMatrix& dm);


string NeighborJoining(DistanceMatrix& dm);


string UPGMA(TaxonSet& ts, DistanceMatrix& dm);
#endif
