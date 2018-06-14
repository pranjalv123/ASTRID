#include "catch.hpp"
#include <newick.hpp>
#include "../octal.hpp"


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
