#include <TreeClade.hpp>
#include <newick.hpp>
#include "octal.hpp"




TreeClade& find_root(TreeClade& T, Taxon x,  const Clade& t_taxa) {
  for (int i = 0; i < T.nchildren(); i++) {
    TreeClade& c = T.child(i);
    if (!c.contains(x))
      continue;
    if (c.overlap_size(t_taxa) == 0)
      return T;
    return find_root(c, x, t_taxa);
  }
}


void insert_leaf(Tree& T, Tree& t, int T_node, int t_node, Taxon x, unordered_map<int, int>& A) {
  int newnode = t.addNode();

  int xnode = t.addNode();

  t.node(newnode) += t.node(t_node);
  t.node(xnode) += x;

  t.node(t_node) += x;


  for (int i = 0; i < t.node(t_node).nchildren(); i++) {
    t.node(newnode).addChild(t.node(t_node).children().at(i));
  }

  t.node(t_node).children().clear();

  t.node(t_node).addChild(newnode);
  t.node(t_node).addChild(xnode);

  int current = t_node;
  while(current >= 0) {
    t.node(current) += x;
    current = t.node(current).parent;
  }

  vector<int> to_change;
  int T_xnode = -1;
  for (auto x : A) {
    if (x.second == t_node && x.first != T_node ){
      to_change.push_back(x.first);
    }
  }
  for (int i : to_change) {
    A[i] = newnode;
  }

  for (int i = 0; i < T.next_entry; i++) {
    if (T.node(i).size() == 1 && T.node(i).contains(x)){
      A[i] = xnode;
      break;
    }
  }



}

void add_node(Tree& T, Tree& t, unordered_map<int, int>& A, Taxon x) {
  int xroot = find_root(T.root(), x, t.taxa()).index;
  while(true) {
    if(A.count(xroot)) {
      insert_leaf(T, t, xroot,  A[xroot], x, A);
      return;
    }
    for (int i = 0; i < T.node(xroot).nchildren(); i++) {
      if (T.node(xroot).child(i).overlap_size(t.taxa()) > 0) {
        xroot = T.node(xroot).child(i).index;
        break;
      }
    }
  }
}

void shared_edges(Tree& T, Tree& t, unordered_map<int, int>& A) {
  TaxonSet& ts = T.ts;
  DistanceMatrix lca1(ts);
  DistanceMatrix lca2(ts);
  T.LCA(lca1);
  t.LCA(lca2);


  for (int i = 0; i < T.clades.size(); i++) {
    TreeClade& c = T.node(i);

    if (c.overlap_size(t.root()) == 0) {
      continue;
    }

    int current = -1;
    int tax1 = -1;
    for (Taxon tax : c) {
      if (!t.root().contains(tax)) {
        continue;
      }
      if (current == -1) {
        current = lca2(tax, tax);
        tax1 = tax;
      }
      else {
        int test = lca2(tax1, tax);
        if (t.node(current).size() < t.node(test).size())
          current = test;
      }

    }

    int matching = 1;
    for (Taxon tax: c) {
      if (!t.root().contains(tax)) {
        continue;
      }
      if (!t.node(current).contains(tax)){
        matching=0;
        break;
      }
    }
    if (matching) {
      A[i] = current;
    } else {
    }
  }
}

void octal_complete(Tree& T, Tree& t) {
  unordered_map<int, int> shared;
  TaxonSet& ts = T.ts;
  shared_edges(T, t, shared);
  for (Taxon x : T.taxa()) {
    if (!t.taxa().contains(x)) {
      add_node(T, t, shared, x);
    }
  }

}
