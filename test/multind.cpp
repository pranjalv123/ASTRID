#include "../src/multind.hpp"
#include "catch2.hpp"
#include "phylokit/TaxonSet.hpp"
#include "phylokit/newick.hpp"
#include <iostream>

TEST_CASE("test averaging") {
  TaxonSet indiv(6);
  indiv.add("A1");
  indiv.add("A2");
  indiv.add("B1");
  indiv.add("C1");
  indiv.add("C2");
  indiv.add("C3");

  IndSpeciesMapping mapping(indiv);
  mapping.add_indiv_to_species(indiv["A1"], "A");
  mapping.add_indiv_to_species(indiv["A2"], "A");
  mapping.add_indiv_to_species(indiv["B1"], "B");
  mapping.add_indiv_to_species(indiv["C1"], "C");
  mapping.add_indiv_to_species(indiv["C2"], "C");
  mapping.add_indiv_to_species(indiv["C3"], "C");

  DistanceMatrix indiv_dm1(indiv, "((A1,B1),(A2,C2),C1);");
  DistanceMatrix indiv_dm2(indiv, "(A1,B1,C3);");

  DistanceMatrix species_dm1 = mapping.average(indiv_dm1);
  DistanceMatrix species_dm2 = mapping.average(indiv_dm2);

  species_dm1 += species_dm2;
  for (Taxon i : mapping.species()) {
    for (Taxon j : mapping.species()) {
      if ( i > j && species_dm1.masked(i, j)) {
        species_dm1(i, j) /= species_dm1.masked(i, j);
      }
    }
  }

  REQUIRE(species_dm1(mapping.species()["A"], mapping.species()["A"]) == 0);
  REQUIRE(species_dm1(mapping.species()["A"], mapping.species()["B"]) == 2.5);
  REQUIRE(species_dm1(mapping.species()["A"], mapping.species()["C"]) == 2.5);
  REQUIRE(species_dm1(mapping.species()["B"], mapping.species()["B"]) == 0);
  REQUIRE(species_dm1(mapping.species()["B"], mapping.species()["C"]) == 2.75);
  REQUIRE(species_dm1(mapping.species()["C"], mapping.species()["C"]) == 0);
}