#ifndef ASTRID_MULTIND_ONCE__
#define ASTRID_MULTIND_ONCE__

#include <fstream>
#include <iostream>
#include <DistanceMatrix.hpp>
#include <string>
#include <unordered_map>


enum file_format {ASTRAL, ASTRIDM};

class IndSpeciesMapping {
private:
  map<Taxon, Taxon> ind_species_map;
  map<Taxon, vector<Taxon>> species_ind_map;

  TaxonSet& indiv_ts;
  TaxonSet species_ts;
  
public:
  IndSpeciesMapping(TaxonSet& indiv_ts_) :
    indiv_ts(indiv_ts_),
    species_ts(indiv_ts.size())
  {}
  
  //IndSpeciesMapping(istream& instream);
  //IndSpeciesMapping(string& infile);

  Taxon operator[](Taxon t);

  DistanceMatrix average(DistanceMatrix& indiv_mat) const;

  file_format identify(istream& instream);

  void load_astral(istream& instream);
  void load_astridm(istream& instream);
  void load(istream& instream);
  void load(string& infile);      
  
  TaxonSet& species();
  TaxonSet& indivs();
  
};

#endif
