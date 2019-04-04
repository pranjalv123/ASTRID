#include "catch.hpp"
#include <newick.hpp>
#include "../octal.hpp"
#include "../Args.hpp"


TEST_CASE("OCTAL") {
  TaxonSet ts(10);
  for (int i = 0; i < 10; i++) {
    stringstream ss;
    ss << i;
    ts.add(ss.str());
  }

  string t1_string = "((1,2),(4,(3,5)))";

  Tree t1 = newick_to_treeclades(t1_string, ts);
  Tree t1_orig = newick_to_treeclades(t1_string, ts);
  Tree t2 = newick_to_treeclades("((1,(2,8)),(3,((6,7),(4,5)))", ts);

  for (int i = 0; i < t1.next_entry; i++) {
    cout << i << ": " << t1.node(i) << endl;
  }

  octal_complete(t2, t1);

  cout << t1 << endl;
  cout << t2 << endl;
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

  REQUIRE(!args.octal);
  
}
