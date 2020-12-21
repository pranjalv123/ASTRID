#include "multind.hpp"
#include <boost/tokenizer.hpp>
#include <glog/logging.h>


file_format IndSpeciesMapping::identify(std::istream& instream) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" :,");
  std::string line;

  file_format ff = ASTRIDM;
  
  while(std::getline(instream, line)) {
    if (line.size() == 0)
      continue;
    if (line.find(":") != std::string::npos) {
      ff = ASTRAL;
      break;
    }
    if (line.find(",") != std::string::npos) {
      ff = ASTRAL;
      break;
    }


    tokenizer tizer(line, sep);
    std::vector<std::string> tokens(tizer.begin(), tizer.end());
    
    if (tokens.size() > 2) {
      ff = ASTRAL;
      break;
    }

    if (indiv_ts.has(tokens[0]) && (! indiv_ts.has(tokens[1]))) {
      ff = ASTRIDM;
      break;
    }

    if (indiv_ts.has(tokens[1]) && (! indiv_ts.has(tokens[0]))) {
      ff = ASTRAL;
      break;
    } 
      
  }
  
  instream.clear();
  instream.seekg(0, std::ios::beg);
  return ff;
}

void IndSpeciesMapping::load_astral(std::istream& instream) {
  VLOG(1) << "Loading ASTRAL mapping file";
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" :,");
  std::string line;
  int lineno = 0;
  
  while (std::getline(instream,line)) {
    lineno++;
    if (line.size() == 0)
      continue;
    tokenizer tokens(line, sep);
    std::vector<std::string> temp;
    temp.assign(tokens.begin(),tokens.end());

    if (temp.size() < 2) {
      LOG(ERROR) << temp.size() << " tokens on line " << lineno << " in mapping file (expected at least 2)" << std::endl;
      LOG(ERROR) << "Line " << lineno << ": " <<  line << std::endl;
      exit(0);
    }
    
    Taxon species_tax = species_ts.add(temp[0]);

    for (size_t i = 1; i < temp.size(); i++) {
      
      if( ! indiv_ts.has(temp[i]) ) {
	      LOG(ERROR) << "Unrecognized taxon " << temp[i] << " in mapping file" << std::endl;
	      LOG(ERROR) << "Line " << lineno << ": " <<  line << std::endl;
	      exit(0);
      }
      
      Taxon indiv_tax   = indiv_ts[temp[i]];    
      
      if (ind_species_map.count(indiv_tax)) {
	      LOG(ERROR) << "Individual " << temp[i] << " appears twice in mapping file!" << std::endl;
	      LOG(ERROR) << "Line " << lineno << ": " <<  line << std::endl;
	      exit(0);
      }

    
      ind_species_map[indiv_tax] = species_tax;
      species_ind_map[species_tax].push_back(indiv_tax);    
      
    }
  }
  VLOG(1) << "Loaded ASTRAL mapping file";
}

void IndSpeciesMapping::load_astridm(std::istream& instream) {
  VLOG(1) << "Loading ASTRIDM mapping file";
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" ");
  std::string line;
  int lineno = 0;
  
  while (std::getline(instream,line)) {
    lineno++;
    if (line.size() == 0)
      continue;
    tokenizer tokens(line, sep);
    
    std::vector<std::string> temp;
    temp.assign(tokens.begin(),tokens.end());
    if (temp.size() != 2) {
      LOG(ERROR) << temp.size() << " tokens on line " << lineno << " in mapping file (expected 2)" << std::endl;
      LOG(ERROR) << "Line " << lineno << ": " <<  line << std::endl;
      exit(0);
    }

    if( ! indiv_ts.has(temp[0]) ) {
      continue; // DO NOT SUBMIT
      LOG(ERROR) << "Unrecognized taxon " << temp[0] << " in mapping file" << std::endl;
      LOG(ERROR) << "Line " << lineno << ": " <<  line << std::endl;
      exit(0);
    }

    
    Taxon species_tax = species_ts.add(temp[1]);
    Taxon indiv_tax   = indiv_ts[temp[0]];    

    if (ind_species_map.count(indiv_tax)) {
      LOG(ERROR) << "Individual " << temp[0] << " appears twice in mapping file!" << std::endl;
      LOG(ERROR) << "Line " << lineno << ": " <<  line << std::endl;
      exit(0);
    }

    
    ind_species_map[indiv_tax] = species_tax;
    species_ind_map[species_tax].push_back(indiv_tax);    

  }
  VLOG(1) << "Loaded ASTRIDM mapping file";
}


void IndSpeciesMapping::load(std::istream& instream) {
  if (identify(instream) == ASTRAL) {
    load_astral(instream);
  } else {
    load_astridm(instream);
  }

  VLOG(1) << "ind_species_map has " << ind_species_map.size() << " members";
  VLOG(1) << "species_ind_map has " << species_ind_map.size() << " members";

}



void IndSpeciesMapping::load(std::string& infile) {
  std::ifstream instream(infile);
  load(instream);
  instream.close();
}

Taxon IndSpeciesMapping::operator[](Taxon t) {
  return ind_species_map.at(t);
}

DistanceMatrix IndSpeciesMapping::average(DistanceMatrix& indiv_mat) const {
  DistanceMatrix species_mat(species_ts);
  for (const Taxon i_s : species_ts) {
    for (const Taxon j_s : species_ts) {
      if (i_s == j_s) {
        species_mat(i_s, j_s) = 0;
        species_mat.masked(i_s, j_s) = 1;
        continue;
      }
      double sum_ij = 0;
      double count_ij = 0;
      for (const Taxon i_i : species_ind_map.at(i_s) ) {
        for (const Taxon j_i : species_ind_map.at(j_s) ) {
          if (indiv_mat.has(i_i, j_i)) {
            sum_ij += indiv_mat(i_i, j_i);
            count_ij++;
          }
        }
      }
      if (count_ij == 0) {
	      species_mat(i_s, j_s) = 0;
	      species_mat.masked(i_s, j_s) = 0;	
      } else {      
	      species_mat(i_s, j_s) = sum_ij/count_ij;
	      species_mat.masked(i_s, j_s) = 1;		
      }
    }
  }
  return species_mat;
}

DistanceMatrix IndSpeciesMapping::average(SparseDistanceMatrix& indiv_mat) const {
  DistanceMatrix species_mat(species_ts);


  for (const Taxon i : species_ts) {
    for (const Taxon j : species_ts) {
      species_mat(i, j) = 0;
      species_mat.masked(i, j) = 0;
    }
  }

  for (auto& entry : indiv_mat) {
    Taxon i_indiv = entry.first.first;
    Taxon j_indiv = entry.first.second;

    Taxon i_species = ind_species_map.at(i_indiv);
    Taxon j_species = ind_species_map.at(j_indiv);

    double distance = entry.second;

    if(i_species == j_species) {
      species_mat(i_species, j_species) = 0;
      species_mat.masked(i_species, j_species) = 1;
      continue;
    }

    species_mat(i_species, j_species) += distance;
    species_mat.masked(i_species, j_species) += 1;
  }

  for (const Taxon i : species_ts) {
    for (const Taxon j : species_ts) {
      
      if (species_mat.masked(i, j) > 0) {
        species_mat(i, j) /= species_mat.masked(i, j);
        species_mat.masked(i, j) = 1;
      }
    }
  }
  return species_mat;
}

void IndSpeciesMapping::add_indiv_to_species(Taxon indiv, std::string species) {
  ind_species_map[indiv] = species_ts.add(species);
  species_ind_map[species_ts.add(species)].push_back(indiv);
}

TaxonSet& IndSpeciesMapping::species() {
  return species_ts;
}
TaxonSet& IndSpeciesMapping::indivs() {
  return indiv_ts;
}
