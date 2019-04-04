#include "DistanceMethods.hpp"
#include "Args.hpp"
#include "octal.hpp"
#include <newick.hpp>
#include <iostream>
#include <fstream>
#include "util/Logger.hpp"


bool has_missing(TaxonSet& ts, DistanceMatrix& dm) {
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i+1; j < ts.size(); j++) {
      if ( ! dm.has(i, j) ) {
	return true;
      }
    }
  }
  return false;
}

TaxonSet get_ts(vector<string>& newicks){
  unordered_set<string> taxa;
  for (string n : newicks) {
    newick_to_ts(n, taxa);
  }
  TaxonSet ts(taxa.size());
  for (string t : taxa) {
    ts.add(t);
  }
  return ts;

 }
DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks, vector<double> weights, vector<Clade>& tree_taxa) {

  DistanceMatrix result(ts);
  for (size_t i = 0; i < newicks.size(); i++) {
    string& n = newicks[i];
    double w  = weights[i];


    DistanceMatrix dm(ts, n);
    if (tree_taxa.size() > i) {
      for (Taxon t1 : ts) {
        for (Taxon t2 : ts) {
          if ( ! (tree_taxa[i].contains(t1) && tree_taxa[i].contains(t2))) {
            dm(t1, t2) = 0;
            dm.masked(t1, t2) = 0;
          }
        }
      }
    }

    dm *= w;
    result += dm;

  }

  for (int i = 0; i < ts.size(); i++) {
    for (int j = i; j < ts.size(); j++) {
        //cout << result(i,j) << "," << result.masked(i,j) << "\t";
        if (result.masked(i,j))
          result(i,j) /= result.masked(i,j);
    }
  }
//  cout << endl;

  return result;

}

DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks) {
  vector<Clade> vc;
  return get_distance_matrix(ts, newicks, vector<double>(newicks.size(), 1), vc);
}

DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks, vector<Clade>& tree_taxa) {
  return get_distance_matrix(ts, newicks, vector<double>(newicks.size(), 1), tree_taxa);
}

string run_astrid(vector<string> newicks) {
  return "";
}


void fill_in(TaxonSet& ts, DistanceMatrix& dm, double cval) {
  int count = 0;
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        dm(i,j) = cval;
        dm.masked(i,j) = 1;
        count++;
      }
    }
  }
  cerr << "Filled in " << count << " elements" << endl;

}

void fill_in(TaxonSet& ts, DistanceMatrix& dm, string tree) {
  vector<string> trees;
  trees.push_back(tree);
  DistanceMatrix dm_tree = get_distance_matrix(ts, trees);

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        dm(i,j) = dm_tree(i,j);
        dm.masked(i,j) = 1;
      }
    }
  }

}

void progressbar(double pct) {
  cerr << "[";
  for (int i = 0; i < (int)(68 * pct); i++) {
    cerr << "#";
  }
  for (int i = (int)(68 * pct); i < 68; i++) {
    cerr << "-";
  }
  cerr << "]\r";
}

int main(int argc, char** argv) {
  Logger::enable("PROGRESS");
  Args args(argc, argv);
  vector<string> input_trees;

  ifstream inf(args.infile);

  string buf;
  PROGRESS << "Reading trees..." << endl;
  while(!inf.eof()) {
    getline(inf, buf);
    if (buf.size() > 3)
      input_trees.push_back(buf);
  }
  PROGRESS << "Read " << input_trees.size() << " trees" << endl;
  
  TaxonSet ts = get_ts(input_trees);
  int iter = 1;
  string tree;

  
  DistanceMatrix dm = get_distance_matrix(ts, input_trees);

  cerr << "Estimating tree" << endl;
  for (string method : args.dms) {
    cerr << "Running " << method << endl;

    if (tree.size() && args.octal) {

      vector<string> completed_trees;
      vector<Clade>  tree_taxa;
      Tree T = newick_to_treeclades(tree, ts);
      for (string t_s : input_trees) {
        Tree t = newick_to_treeclades(t_s, ts);
        tree_taxa.push_back(t.taxa());

        octal_complete(T, t);
        stringstream ss;
        ss << t;
        completed_trees.push_back(ss.str());
      }
      dm = get_distance_matrix(ts, input_trees, tree_taxa);

    }

    if (iter > 1) {
      fill_in(ts, dm, tree);
    } else if (args.constant != 0) {
      fill_in(ts, dm, args.constant);
    }

    if (method == "auto") {
      if (has_missing(ts, dm)) {
	tree = BioNJStar(ts, dm);
      } else {
	tree = FastME(ts, dm, 1, 1);
      }
    } else if (method == "upgma") {
      tree = UPGMA(ts, dm);
    } else if (method == "fastme") {
      tree = FastME(ts, dm, 0, 0);
    } else if (method == "fastme_nni") {
      tree = FastME(ts, dm, 1, 0);
    } else if (method == "fastme_spr") {
      tree = FastME(ts, dm, 1, 1);
    } else if (method == "bionj") {
      tree = BioNJStar(ts, dm);
    }

    ofstream outfile(args.outfile + "." + to_string(iter));
    outfile << tree << endl;
    iter ++;
  }
  ofstream outfile(args.outfile);
  outfile << tree << endl;

  cout <<  tree << endl;

  return 0;
}
