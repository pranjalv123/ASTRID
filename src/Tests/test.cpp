#include "catch.hpp"
#include <newick.hpp>
#include "../octal.hpp"
#include "../expand.hpp"
#include <iostream>

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



TEST_CASE("remove_taxa") {
  TaxonSet ts(10);
  for (int i = 0; i < 10; i++) {
    stringstream ss;
    ss << i;
    ts.add(ss.str());
  }

  string t1_string = "(8,9,((1,2),(4,(3,(5,6,7)))))";

  Clade tree_taxa(ts, ts.taxa_bs);

  tree_taxa.remove(9);
  tree_taxa.remove(6);

  tree_taxa.remove(1);
//  tree_taxa.remove(2);
  string rmt = remove_missing_taxa(ts, tree_taxa, t1_string);

  cout << t1_string << endl;

  cout << rmt << endl;
}

TEST_CASE("expand") {
  TaxonSet ts(10);
  for (int i = 0; i < 10; i++) {
    stringstream ss;
    ss << i;
    ts.add(ss.str());
  }

  string t1_string = "((1,2), ((3,6), (4, 5)))";
  string t2_string = "((1, 2), (3,4)";

  DistanceMatrix dm(ts, t2_string);
  DistanceMatrix dm_guide(ts, t1_string);

  for (Taxon t1 : ts) {
    for (Taxon t2 : ts) {
      if (t1 > 0 && t1 < 5 && t2 > 0 && t2 < 5)
        cout << dm(t1, t2) << " ";
    }
    cout << endl;
  }
  cout << endl << endl;

  for (Taxon t1 : ts) {
    for (Taxon t2 : ts) {
      if (t1 > 0 && t1 < 5 && t2 > 0 && t2 < 5)
        cout << dm_guide(t1, t2) << " ";
    }
    cout << endl;
  }
  cout << endl << endl;


  Clade tt(ts);
  tt.add(1);
  tt.add(2);
  tt.add(3);
  tt.add(4);

  do_expand(dm, ts,tt, t2_string, dm_guide);


  for (Taxon t1 : ts) {
    for (Taxon t2 : ts) {
      if (t1 > 0 && t1 < 5 && t2 > 0 && t2 < 5)
        cout << dm(t1, t2) << " ";
    }
    cout << endl;
  }





}
