#include "DistanceMethods.hpp"
#include "Args.hpp"
#include "octal.hpp"
#include "expand.hpp"
#include <newick.hpp>
#include <iostream>
#include <fstream>
#include <sstream>

TaxonSet get_ts(vector<string> newicks){
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

DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks, vector<double> weights, vector<Clade>& tree_taxa, bool expand, string guide_tree) {
  DistanceMatrix guide_dm(ts);
  if (expand && guide_tree.size()) {
      guide_dm = DistanceMatrix(ts, guide_tree);
      cout << "GOT GUIDE DM" << endl;
    }

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


    if (expand) {
      do_expand(dm, ts, tree_taxa[i], guide_tree, guide_dm);
    }


    result += dm;

  }

  for (int i = 0; i < ts.size(); i++) {
    for (int j = i; j < ts.size(); j++) {
        if (result.masked(i,j))
          result(i,j) /= result.masked(i,j);
    }
  }

  return result;

}

DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks) {
  vector<Clade> vc;
  return get_distance_matrix(ts, newicks, vector<double>(newicks.size(), 1), vc, false, "");
}

DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks, vector<Clade>& tree_taxa, bool expand, string guide_tree) {
  return get_distance_matrix(ts, newicks, vector<double>(newicks.size(), 1), tree_taxa, expand, guide_tree);
}


DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks, vector<Clade>& tree_taxa) {
  return get_distance_matrix(ts, newicks, vector<double>(newicks.size(), 1), tree_taxa, false, "");
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

int main(int argc, char** argv) {
  Args args(argc, argv);
  vector<string> input_trees;

  ifstream inf(args.infile);

  string buf;
  while(!inf.eof()) {
    getline(inf, buf);
    if (buf.size() > 3)
      input_trees.push_back(buf);
  }

  TaxonSet ts = get_ts(input_trees);
  int iter = 1;
  string tree;
  DistanceMatrix dm = get_distance_matrix(ts, input_trees);


  for (string method : args.dms) {
    cout << "Running " << method << endl;

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

    if (tree.size() && args.expand) {

      vector<string> completed_trees;
      vector<Clade>  tree_taxa;
      for (string t_s : input_trees) {
        Tree t = newick_to_treeclades(t_s, ts);
        tree_taxa.push_back(t.taxa());
      }
      dm = get_distance_matrix(ts, input_trees, tree_taxa, true, tree);
    }

    if (iter > 1) {
      fill_in(ts, dm, tree);
    } else if (args.constant != 0) {
      fill_in(ts, dm, args.constant);
    }


    if (method == "upgma") {
      tree = UPGMA(ts, dm);
    } else if (method == "fastme") {
      tree = FastME(ts, dm, 0, 0);
    } else if (method == "fastme_nni") {
      tree = FastME(ts, dm, 1, 0);
    } else if (method == "fastme_spr") {
      tree = FastME(ts, dm, 1, 1);
    }

    ofstream outfile(args.outfile + "." + to_string(iter));
    outfile << tree << endl;
    iter ++;
  }
  ofstream outfile(args.outfile);
  outfile << tree << endl;

  cout <<  tree << endl;


}
