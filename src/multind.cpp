#include "multind.hpp"
#include <util/Logger.hpp>
#include<boost/tokenizer.hpp>


file_format IndSpeciesMapping::identify(istream& instream) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" :,");
  string line;

  file_format ff = ASTRIDM;
  
  while(getline(instream, line)) {
    if (line.size() == 0)
      continue;
    if (line.find(":") != string::npos) {
      ff = ASTRAL;
      break;
    }
    if (line.find(",") != string::npos) {
      ff = ASTRAL;
      break;
    }


    tokenizer tizer(line, sep);
    vector<string> tokens(tizer.begin(), tizer.end());
    
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
  instream.seekg(0, ios::beg);
  return ff;
}

void IndSpeciesMapping::load_astral(istream& instream) {

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" :,");
  string line;
  int lineno = 0;
  
  while (getline(instream,line)) {
    lineno++;
    if (line.size() == 0)
      continue;
    tokenizer tokens(line, sep);
    vector<string> temp;
    temp.assign(tokens.begin(),tokens.end());

    if (temp.size() < 2) {
      ERR << temp.size() << " tokens on line " << lineno << " in mapping file (expected at least 2)" << endl;
      ERR << "Line " << lineno << ": " <<  line << endl;
      exit(0);
    }
    
    Taxon species_tax = species_ts.add(temp[0]);

    for (int i = 1; i < temp.size(); i++) {
      
      if( ! indiv_ts.has(temp[i]) ) {
	ERR << "Unrecognized taxon " << temp[i] << " in mapping file" << endl;
	ERR << "Line " << lineno << ": " <<  line << endl;	
	exit(0);
      }
      
      Taxon indiv_tax   = indiv_ts[temp[i]];    
      
      if (ind_species_map.count(indiv_tax)) {
	ERR << "Individual " << temp[i] << " appears twice in mapping file!" << endl;
	ERR << "Line " << lineno << ": " <<  line << endl;
	exit(0);
      }

    
      ind_species_map[indiv_tax] = species_tax;
      species_ind_map[species_tax].push_back(indiv_tax);    
      
    }
  }
}

void IndSpeciesMapping::load_astridm(istream& instream) {

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" ");
  string line;
  int lineno = 0;
  
  while (getline(instream,line)) {
    lineno++;
    if (line.size() == 0)
      continue;
    tokenizer tokens(line, sep);
    
    vector<string> temp;
    temp.assign(tokens.begin(),tokens.end());
    if (temp.size() != 2) {
      ERR << temp.size() << " tokens on line " << lineno << " in mapping file (expected 2)" << endl;
      ERR << "Line " << lineno << ": " <<  line << endl;
      exit(0);
    }

    if( ! indiv_ts.has(temp[0]) ) {
      ERR << "Unrecognized taxon " << temp[0] << " in mapping file" << endl;
      ERR << "Line " << lineno << ": " <<  line << endl;	
      exit(0);
    }

    
    Taxon species_tax = species_ts.add(temp[1]);
    Taxon indiv_tax   = indiv_ts[temp[0]];    

    if (ind_species_map.count(indiv_tax)) {
      ERR << "Individual " << temp[0] << " appears twice in mapping file!" << endl;
      ERR << "Line " << lineno << ": " <<  line << endl; 
      exit(0);
    }

    
    ind_species_map[indiv_tax] = species_tax;
    species_ind_map[species_tax].push_back(indiv_tax);    

  }
  return;
}


void IndSpeciesMapping::load(istream& instream) {
  if (identify(instream) == ASTRAL) {
    load_astral(instream);
  } else {
    load_astridm(instream);
  }

  for (Taxon indiv_tax : indiv_ts) {
    if (ind_species_map.count(indiv_tax) == 0) {
      Taxon species_tax = species().add(indiv_ts[indiv_tax]);

      ind_species_map[indiv_tax] = species_tax;
      species_ind_map[species_tax].push_back(indiv_tax);
    }
  }
}



void IndSpeciesMapping::load(string& infile) {
  ifstream instream(infile);
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
      double sum_ij = 0;
      double count_ij = 0;
      for (const Taxon i_i : species_ind_map.at(i_s) ) {
        for (const Taxon j_i : species_ind_map.at(j_s) ) {
	    sum_ij += indiv_mat(i_i, j_i);
	    count_ij++;
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

TaxonSet& IndSpeciesMapping::species() {
  return species_ts;
}
TaxonSet& IndSpeciesMapping::indivs() {
  return indiv_ts;
}
