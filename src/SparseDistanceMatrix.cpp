#include "SparseDistanceMatrix.hpp"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <glog/logging.h>


using namespace std;

SparseDistanceMatrix::SparseDistanceMatrix(const TaxonSet& ts) : ts(ts){
}

/**
* Add all the members from d with non-zero mask
**/
SparseDistanceMatrix::SparseDistanceMatrix(const TaxonSet& ts, DistanceMatrix& d) 
: ts(ts) {
    for (Taxon i : ts) {
        for(Taxon j : ts) {
            if (d.has(i, j)) {
                set(i, j, d.get(i,j));
            }
        }
    }
}

SparseDistanceMatrix::SparseDistanceMatrix(const TaxonSet& ts, const string& newick) : ts(ts)  {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  tokenizer tokens(newick, sep);

  std::vector<double> dists(ts.size(), 0);
  std::vector<double> distance_to_lca(ts.size(), 0);

  std::vector<Taxon> seen(ts.size(), 0);
  std::string prevtok = "";

  for (auto tok : tokens) {
    VLOG_EVERY_N(1, 100) << "Processed " << google::COUNTER << " tokens";
    if (tok == "(") {
      for (Taxon s : seen) {
        distance_to_lca[s] += 1;
        distance_to_lca[s] += 1;
      }
    } else if (tok == ")") {
      for (Taxon s : seen) {
        if (distance_to_lca[s] == 0) { // we're in a sister leaf to s
          dists[s] += 1;
        } else { 
          dists[s] -= 1;
          distance_to_lca[s] -= 1;
        } 
      }
    } else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if (prevtok == ")" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        continue;
      }
      VLOG_EVERY_N(1, 100) << "Processed " << google::COUNTER << " taxa";

      boost::algorithm::trim(tok);
      Taxon id = ts[tok];
      for (Taxon other : seen) {
        if (contains(other, id)) {
            set(other, id, get(other, id) + dists[other] + 2);
        } else {
            set(other, id, dists[other] + 2);
        }
      }
      seen.push_back(id);
    }
    prevtok = tok;
  }
}

/** 
* Set the distance between two taxa
*/
void SparseDistanceMatrix::set(Taxon t1, Taxon t2, double distance) {
    Taxon a = std::min(t1, t2);
    Taxon b = std::max(t1, t2);

    d[make_pair(a, b)] = distance;
    
}

double SparseDistanceMatrix::get(Taxon t1, Taxon t2) const {
    Taxon a = std::min(t1, t2);
    Taxon b = std::max(t1, t2);
    return d.at(make_pair(a, b));
}

bool SparseDistanceMatrix::contains(Taxon t1, Taxon t2) const {
    Taxon a = std::min(t1, t2);
    Taxon b = std::max(t1, t2);
    return d.count(make_pair(a,b)) > 0;
}