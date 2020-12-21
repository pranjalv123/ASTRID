#include "Args.hpp"
#include "DistanceMethods/DistanceMethods.hpp"
#include "multind.hpp"
#include "octal.hpp"
#include "SparseDistanceMatrix.hpp"
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
    VLOG_EVERY_N(1, 100) << "Got taxa from " << google::COUNTER << " trees";
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
                                   IndSpeciesMapping *imap,
                                   bool sparse) {

  DistanceMatrix result = imap ? DistanceMatrix(imap->species()) : DistanceMatrix(result);

  for (size_t i = 0; i < newicks.size(); i++) {
    LOG(INFO) << "derooting";
    std::string newick_derooted = deroot(newicks[i]);
    LOG(INFO) << "derooted";
    double w = weights[i];
    if (sparse && imap) {
      // use a sparse distance matrix for the individuals
      LOG(INFO) << "getting sparse dm";
      SparseDistanceMatrix sdm(ts, newick_derooted);
      LOG(INFO) << "got it; averaging";
      result += imap->average(sdm);
      LOG(INFO) << "averaged";
    } else{
      LOG(INFO) << "getting dm";
      DistanceMatrix dm(ts, newick_derooted);
      LOG(INFO) << "got it";
      if (w != 1) {
        LOG(INFO) << "weighting";
        dm *= w;
      }
      if (imap) {
        LOG(INFO) << "averaging";
        dm = imap->average(dm);
      }
      LOG(INFO) << "adding";
      result += dm;
    }

    VLOG_EVERY_N(1, 1) << "Got " << google::COUNTER << " distance matrices";
    
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
                                   IndSpeciesMapping *imap,
                                   bool sparse) {
  std::vector<Clade> vc;
  return get_distance_matrix(ts, newicks,
                             std::vector<double>(newicks.size(), 1), imap, sparse);
}

std::string run_astrid(std::vector<std::string> newicks) { return ""; }

void fill_in(TaxonSet &ts, DistanceMatrix &dm, double cval, bool sparse) {
  int count = 0;
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        dm(i, j) = cval;
        dm.masked(i, j) = 1;
        count++;
      }
    }
  }
  std::cerr << "Filled in " << count << " elements" << std::endl;
}

void fill_in(TaxonSet &ts, DistanceMatrix &dm, std::string tree, bool sparse) {
  std::vector<std::string> trees;
  trees.push_back(tree);
  DistanceMatrix dm_tree = get_distance_matrix(ts, trees, NULL, sparse);

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        dm(i, j) = dm_tree(i, j);
        dm.masked(i, j) = 1;
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
  FLAGS_logtostderr = 1;
  google::InitGoogleLogging(argv[0]);

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

  LOG(INFO) << "Getting taxa from trees...";
  TaxonSet ts = get_ts(input_trees);
  LOG(INFO) << "Found " << ts.size() << " taxa" << std::endl;

  // Set up multi-individual mapping
  IndSpeciesMapping *multind_mapping = nullptr;
  if (args.multindfile.size()) {
    multind_mapping = new IndSpeciesMapping(ts);
    multind_mapping->load(args.multindfile);
  }

  int iter = 1;
  std::string tree;

  DistanceMatrix dm = get_distance_matrix(ts, input_trees, multind_mapping, args.sparse_individual_matrix);

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
      dm = get_distance_matrix(ts, input_trees, multind_mapping, args.sparse_individual_matrix);
    }

    TaxonSet *species_ts = &ts;
    if (multind_mapping) {
      species_ts = &(multind_mapping->species());
    }

    // fill in missing elements on second and later iterations
    if (iter > 1) {
      fill_in(*species_ts, dm, tree, args.sparse_individual_matrix);
    } else if (args.constant != 0) {
      fill_in(*species_ts, dm, args.constant, args.sparse_individual_matrix);
    }

    if (method == "auto") {
      if (has_missing(*species_ts, dm)) {
        std::cerr << "Missing entries in distance matrix, trying to run BioNJ*"
                  << std::endl;
        std::cerr << "You may have better results adding -u -s to your command "
                     "line to use UPGMA completion instead."
                  << std::endl;
        tree = BioNJStar(*species_ts, dm, args.java_opts);
      } else {
        std::cerr
            << "No missing entries in distance matrix, running FastME2+SPR"
            << std::endl;
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

    std::ofstream outfile(args.outfile + "." + std::to_string(iter));
    outfile << tree << std::endl;
    iter++;
  }
  std::ofstream outfile(args.outfile);
  outfile << tree << std::endl;

  std::cout << tree << std::endl;

  return 0;
}
