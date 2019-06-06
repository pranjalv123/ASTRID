#ifndef ASTRID_DISTANCE_METHODS__
#define ASTRID_DISTANCE_METHODS__


#include <DistanceMatrix.hpp>
#include <queue>
#include <sstream>


string FastME (TaxonSet& ts, DistanceMatrix& dm, int nni, int spr);

string RapidNJ (TaxonSet& ts, DistanceMatrix& dm);

string BioNJStar(TaxonSet& ts, DistanceMatrix& dm, vector<string>& java_opts);


string NeighborJoining(DistanceMatrix& dm);


string UPGMA(TaxonSet& ts, DistanceMatrix& dm);
#endif
