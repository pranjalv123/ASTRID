#include "catch2.hpp"
#include "phylokit/newick.hpp"
#include "../src/octal.hpp"
#include "../src/Args.hpp"
#include "../src/multind.hpp"




TEST_CASE("OCTAL") {
  TaxonSet ts(10);
  for (int i = 0; i < 10; i++) {
    std::stringstream ss;
    ss << i;
    ts.add(ss.str());
  }

  std::string t1_string = "((1,2),(4,(3,5)))";

  Tree t1 = newick_to_treeclades(t1_string, ts);
  Tree t1_orig = newick_to_treeclades(t1_string, ts);
  Tree t2 = newick_to_treeclades("((1,(2,8)),(3,((6,7),(4,5)))", ts);


  octal_complete(t2, t1);

  REQUIRE(t1_orig.RFDist(t1) == 0);
  REQUIRE(t1_orig.RFDist(t2, false) == t1.RFDist(t2, false));
}


TEST_CASE("ARGS") {
  char* argv[] = {"ASTRID", "-i", "inputfile", "-o", "outputfile", "-u", "-f", "-n", "-s", "-x", "2.3"};
  Args args(sizeof(argv)/sizeof(*argv), argv);
  
  REQUIRE(args.infile == "inputfile");
  REQUIRE(args.outfile == "outputfile");
  REQUIRE(args.constant == 2.3);
  REQUIRE(args.dms.size() == 4);
  REQUIRE(args.dms[0] == "upgma");
  REQUIRE(args.dms[1] == "fastme");
  REQUIRE(args.dms[2] == "fastme_nni");
  REQUIRE(args.dms[3] == "fastme_spr");
  
}


void verify_mapping(std::string mapfile, TaxonSet& indiv_ts) {
  std::stringstream stream(mapfile);
  IndSpeciesMapping mapping1(indiv_ts);
  mapping1.load(stream);
  REQUIRE(mapping1[indiv_ts["indiv1"]] == mapping1.species()["species1"]);
  REQUIRE(mapping1[indiv_ts["indiv2"]] == mapping1.species()["species1"]);
  REQUIRE(mapping1[indiv_ts["indiv3"]] == mapping1.species()["species2"]);  

}

TEST_CASE("DISCONNECTED") {
  TaxonSet ts(8);
  ts.add("t0");
  ts.add("t1");
  ts.add("t2");
  ts.add("t3");
  ts.add("t4");
  ts.add("t5");
  ts.add("t6");
  ts.add("t7");


std::string t1 = "(t0, t1,(t2,t3));";
std::string t2 = "(t4, t5,(t6,t7));";
  DistanceMatrix dm1(ts, t1);
  DistanceMatrix dm2(ts, t2);

  dm1 += dm2;

  //  cout << UPGMA(ts, dm1);
  
  
}


TEST_CASE("IDENTIFY") {
  TaxonSet indiv_ts(3);
  indiv_ts.add("indiv1");
  indiv_ts.add("indiv2");
  indiv_ts.add("indiv3");

  IndSpeciesMapping mapping(indiv_ts);


std::string astridm_mapfile =
    "indiv1 species1\n"
    "indiv2 species1\n"
    "indiv3 species2\n";
std::stringstream astridm_stream(astridm_mapfile);
  
  REQUIRE(mapping.identify(astridm_stream) == ASTRIDM);
  verify_mapping(astridm_mapfile, indiv_ts);

std::string astral_mapfile_1 =
    "species1 indiv1\n"
    "species1 indiv2\n"
    "species2 indiv3\n";
std::stringstream astral_stream_1(astral_mapfile_1);
  
  REQUIRE(mapping.identify(astral_stream_1) == ASTRAL);
  verify_mapping(astral_mapfile_1, indiv_ts);


std::string astral_mapfile_2 =
    "species1 indiv1 indiv2\n"
    "species2 indiv3\n";
std::stringstream astral_stream_2(astral_mapfile_2);
  
  REQUIRE(mapping.identify(astral_stream_2) == ASTRAL);
  verify_mapping(astral_mapfile_2, indiv_ts);


std::string astral_mapfile_3 =
    "species1:indiv1, indiv2\n"
    "species2:indiv3\n";
std::stringstream astral_stream_3(astral_mapfile_3);
  
  REQUIRE(mapping.identify(astral_stream_3) == ASTRAL);
  verify_mapping(astral_mapfile_3, indiv_ts);
  
  
  
}
