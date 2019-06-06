#include "DistanceMethods.hpp"
#include "multind.hpp"
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

DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks, vector<double> weights, vector<Clade>& tree_taxa, IndSpeciesMapping* imap) {
  
  DistanceMatrix result(ts);
  if (imap) {
    result = DistanceMatrix(imap->species());
  }
  
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

    if (imap) {
      dm = imap->average(dm);
    }
    
    result += dm;

  }

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
        if (result.masked(i,j))
          result(i,j) /= result.masked(i,j);
    }
  }


  return result;

}

DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks, IndSpeciesMapping* imap) {
  vector<Clade> vc;
  return get_distance_matrix(ts, newicks, vector<double>(newicks.size(), 1), vc, imap);
}

DistanceMatrix get_distance_matrix(TaxonSet& ts, vector<string> newicks, vector<Clade>& tree_taxa, IndSpeciesMapping* imap) {
  return get_distance_matrix(ts, newicks, vector<double>(newicks.size(), 1), tree_taxa, imap);
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
  DistanceMatrix dm_tree = get_distance_matrix(ts, trees, NULL);

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
  Logger::enable("ERROR"); 
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
 

  //Set up multi-individual mapping
  IndSpeciesMapping* multind_mapping = nullptr;
  if (args.multindfile.size()) {
    multind_mapping = new IndSpeciesMapping(ts);
    multind_mapping->load(args.multindfile);
  }

  
  DistanceMatrix dm = get_distance_matrix(ts, input_trees, multind_mapping);

  cerr << "Estimating tree" << endl;
  for (string method : args.dms) {
    cerr << "Running " << method << endl;

    //OCTAL completion of trees with missing data
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
      dm = get_distance_matrix(ts, input_trees, multind_mapping);

    }

    TaxonSet* species_ts = &ts;
    if (multind_mapping) {
      species_ts = &(multind_mapping->species());
    } 
    

    //fill in missing elements on second and later iterations
    if (iter > 1) {
      fill_in(*species_ts, dm, tree);
    } else if (args.constant != 0) {
      fill_in(*species_ts, dm, args.constant);
    }

    if (method == "auto") {
      if (has_missing(*species_ts, dm)) {
	cerr << "Missing entries in distance matrix, running BioNJ*" << endl;
	tree = BioNJStar(*species_ts, dm, args.java_opts);
      } else {
	cerr << "No missing entries in distance matrix, running FastME2+SPR" << endl;	
	tree = FastME(*species_ts, dm, 1, 1);
      }
    } else if (method == "upgma") {
      tree = UPGMA(*species_ts, dm);
    } else if (method == "fastme") {
      tree = FastME(*species_ts, dm, 0, 0);
    } else if (method == "fastme_nni") {
      tree = FastME(*species_ts, dm, 1, 0);
    } else if (method == "fastme_spr") {
      tree = FastME(*species_ts, dm, 1, 1);
    } else if (method == "bionj") {
      tree = BioNJStar(*species_ts, dm, args.java_opts);
    } else if (method == "rapidnj") {
      tree = RapidNJ(*species_ts, dm);
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
