#ifndef ASTRID_DISTANCE_METHODS__
#define ASTRID_DISTANCE_METHODS__


#include "phylokit/DistanceMatrix.hpp"
#include <queue>
#include <sstream>


std::string FastME (TaxonSet& ts, DistanceMatrix& dm, int nni, int spr);

std::string RapidNJ (TaxonSet& ts, DistanceMatrix& dm);

std::string BioNJStar(TaxonSet& ts, DistanceMatrix& dm, std::vector<std::string>& java_opts);


std::string NeighborJoining(DistanceMatrix& dm);


std::string UPGMA(TaxonSet& ts, DistanceMatrix& dm);
#endif
