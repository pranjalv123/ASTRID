#include "Args.hpp"
#include "DistanceMethods/DistanceMethods.hpp"
#include "multind.hpp"
#include "octal.hpp"
#include "phylokit/newick.hpp"
#include <fstream>
#include <glog/logging.h>
#include <iostream>

bool has_missing(TaxonSet &ts, DistanceMatrix &dm) {
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i + 1; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        return true;
      }
    }
  }
  return false;
}

TaxonSet get_ts(std::vector<std::string> &newicks) {
  std::unordered_set<std::string> taxa;
  for (std::string n : newicks) {
    newick_to_ts(n, taxa);
  }
  TaxonSet ts(taxa.size());
  for (std::string t : taxa) {
    ts.add(t);
  }
  return ts;
}

DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   std::vector<std::string> newicks,
                                   std::vector<double> weights,
                                   std::vector<Clade> &tree_taxa,
                                   IndSpeciesMapping *imap) {

  DistanceMatrix result(ts);
  if (imap) {
    result = DistanceMatrix(imap->species());
  }

  for (size_t i = 0; i < newicks.size(); i++) {
    std::string newick_derooted = deroot(newicks[i]);
    double w = weights[i];

    DistanceMatrix dm(ts, newick_derooted);
    if (tree_taxa.size() > i) {
      for (Taxon t1 : ts) {
        for (Taxon t2 : ts) {
          if (!(tree_taxa[i].contains(t1) && tree_taxa[i].contains(t2))) {
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
      if (result.masked(i, j))
        result(i, j) /= result.masked(i, j);
    }
  }

  return result;
}

DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   std::vector<std::string> newicks,
                                   IndSpeciesMapping *imap) {
  std::vector<Clade> vc;
  return get_distance_matrix(ts, newicks,
                             std::vector<double>(newicks.size(), 1), vc, imap);
}

DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   std::vector<std::string> newicks,
                                   std::vector<Clade> &tree_taxa,
                                   IndSpeciesMapping *imap) {
  return get_distance_matrix(
      ts, newicks, std::vector<double>(newicks.size(), 1), tree_taxa, imap);
}

void fill_in_const(DistanceMatrix& output, TaxonSet &ts,
                  const DistanceMatrix &dm, double cval) {
  int count = 0;
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        output(i, j) = cval;
        output.masked(i, j) = 1;
        count++;
      }
    }
  }
  std::cerr << "Filled in " << count << " elements" << std::endl;
}

void fill_in(DistanceMatrix &output, TaxonSet &ts, 
             const DistanceMatrix &dm, std::string tree) {
  std::vector<std::string> trees;
  trees.push_back(tree);
  DistanceMatrix dm_tree = get_distance_matrix(ts, trees, NULL);

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        output(i, j) = dm_tree(i, j);
        output.masked(i, j) = 1;
      }
    }
  }
}

void progressbar(double pct) {
  std::cerr << "[";
  for (int i = 0; i < (int)(68 * pct); i++) {
    std::cerr << "#";
  }
  for (int i = (int)(68 * pct); i < 68; i++) {
    std::cerr << "-";
  }
  std::cerr << "]\r";
}

int main(int argc, char **argv) {
  Args args(argc, argv);

  std::vector<std::string> input_trees;
  std::ifstream inf(args.infile);

  std::string buf;
  LOG(INFO) << "Reading trees..." << std::endl;
  while (!inf.eof()) {
    getline(inf, buf);
    if (buf.size() > 3)
      input_trees.push_back(buf);
  }
  LOG(INFO) << "Read " << input_trees.size() << " trees" << std::endl;

  TaxonSet ts = get_ts(input_trees);
  int iter = 1;
  std::string tree;

  // Set up multi-individual mapping
  IndSpeciesMapping *multind_mapping = nullptr;
  if (args.multindfile.size()) {
    multind_mapping = new IndSpeciesMapping(ts);
    multind_mapping->load(args.multindfile);
  }

  DistanceMatrix dm = get_distance_matrix(ts, input_trees, multind_mapping);

  std::cerr << "Estimating tree" << std::endl;
  for (std::string method : args.dms) {
    std::cerr << "Running " << method << std::endl;

    // OCTAL completion of trees with missing data
    if (tree.size() && args.octal) {

      std::vector<std::string> completed_trees;
      std::vector<Clade> tree_taxa;
      Tree T = newick_to_treeclades(tree, ts);
      for (std::string t_s : input_trees) {
        Tree t = newick_to_treeclades(t_s, ts);
        tree_taxa.push_back(t.taxa());

        octal_complete(T, t);
        std::stringstream ss;
        ss << t;
        completed_trees.push_back(ss.str());
      }
      dm = get_distance_matrix(ts, input_trees, multind_mapping);
    }

    TaxonSet *species_ts = &ts;
    if (multind_mapping) {
      species_ts = &(multind_mapping->species());
    }

    bool missing = has_missing(ts, dm);
    DistanceMatrix *filledDM = &dm;
    bool delete_filledDM = false;
    // fill in missing elements on second and later iterations
    if (iter > 1 && missing) {
      filledDM = new DistanceMatrix(ts);
      fill_in(*filledDM, *species_ts, dm, tree);
      delete_filledDM = true;
    } else if (args.constant != 0) {
      filledDM = new DistanceMatrix(ts);
      fill_in_const(*filledDM, *species_ts, dm, args.constant);
      delete_filledDM = true;
    }
    
    if (method == "upgma") {
      tree = UPGMA(*species_ts, *filledDM);
    } else if (method == "fastme") {
      tree = FastME(*species_ts, *filledDM, 0, 0);
    } else if (method == "fastme_nni") {
      tree = FastME(*species_ts, *filledDM, 1, 0);
    } else if (method == "fastme_spr") {
      tree = FastME(*species_ts, *filledDM, 1, 1);
    }

    std::ofstream outfile(args.outfile + "." + std::to_string(iter));
    outfile << tree << std::endl;
    if (args.cache) {
      std::ofstream outfile_cache(args.cachefile + "." + std::to_string(iter));
      filledDM->writePhylip(outfile_cache);
    }
    if (delete_filledDM) {
      delete filledDM;
    }
    iter++;
  }
  std::ofstream outfile(args.outfile);
  outfile << tree << std::endl;

  std::cout << tree << std::endl;

  return 0;
}
