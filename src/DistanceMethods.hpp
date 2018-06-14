#ifndef ASTRID_DISTANCE_METHODS__
#define ASTRID_DISTANCE_METHODS__

extern "C" {
#include "fastme/fastme.h"
}
typedef set mySet;
#include <DistanceMatrix.hpp>
#include <queue>
#include <sstream>


string FastME (TaxonSet& ts, DistanceMatrix& dm, int nni, int spr)
{
	int size = ts.size();
	int numSpecies = size;
	double **D, **A;

	A = initDoubleMatrix (2*numSpecies-2);
	D = initDoubleMatrix (2*numSpecies-2);
	fillZeroMatrix (&A, 2*numSpecies-2);


	for (int i=0; i<size; i++)
	{
		D[i] = (double *) mCalloc (size, sizeof (double));
	}

	for (int i = 0; i < size; i++) {
	  for (int j = 0; j < size; j++) {
	    if (!dm.isMasked(i,j))
	      D[i][j] = dm(i,j);
	  }
	}


	Options options;
	Set_Defaults_Input (&options);
	options.method = TaxAddBAL;
	options.use_SPR = spr;
	options.use_NNI = nni;
	options.NNI    = BALNNI;

	mySet species;
	species.firstNode = 0;
	species.secondNode = 0;

	for (Taxon t : ts) {

		node*  v = makeNode (ts[t].c_str(), -1);
		v->index2 = t;
		addToSet(v, &species);
	}


	tree* t = ComputeTree (&options, D, A, &species, numSpecies, 8);
	int nnicount;
	int sprcount;
	t = ImproveTree(&options, t, D, A, &nnicount, &sprcount, options.fpO_stat_file);

	char* tree_output = new char[size << 10];
	tree_output[0] = '\0';

	NewickPrintTreeStr (t, tree_output, 2);
	return string(tree_output);
}


string BioNJStar(DistanceMatrix& dm) {
  return "";
}


string NeighborJoining(DistanceMatrix& dm) {
  return "";
}

struct Cluster {
	int c1;
	int c2;
	size_t sz;
	unordered_map<size_t, double> d;
	int exists;
	Taxon t;
	Cluster(Taxon t) : c1(-1), c2(-1), sz(1), exists(1), t(t) {}
	Cluster(size_t c1, size_t c2, size_t sz) : c1(c1), c2(c2), sz(sz), exists(1),  t(-1) {}
	priority_queue<pair<double, size_t>, std::vector<pair<double, size_t>>, std::greater<pair<double, size_t>> > pq;

	stringstream newick(TaxonSet& ts, vector<Cluster>& clusts) {
		stringstream ss;
//		cout << c1 << "\t" << c2 << endl;
//		cout << sz << endl;
		if (c1 != -1) {
			ss << "(" << clusts[c1].newick(ts, clusts).str() << "," << clusts[c2].newick(ts, clusts).str() << ")";

		} else {
			ss << ts[t];
		}
		return ss;
	}
};



string UPGMA(TaxonSet& ts, DistanceMatrix& dm) {
	vector<Cluster> clusters;
	for (size_t i = 0; i < ts.size(); i++) {
		clusters.emplace_back(i);
	}

	for (Cluster& c1 : clusters) {
		Taxon i = c1.t;
		for (Cluster& c2 : clusters) {
			Taxon j = c2.t;
			if (i != j && dm.has(i,j)) {
				c1.pq.emplace(dm.get(i, j), j);
				c1.d[j] = dm.get(i,j);
			}
		}
	}

	while(clusters.back().sz < ts.size()) {
		size_t best_cluster = -1;
		double best_score   = numeric_limits<double>::max();

		for (size_t i = 0; i < clusters.size(); i++) {
			Cluster& c = clusters[i];
			if(c.exists) {

				while (!c.pq.empty() && !clusters[c.pq.top().second].exists) {
						//cout << c.pq.top().second << "\t";
					c.pq.pop();
				}
				if (c.pq.empty()) {

					continue;
				}

				if (c.pq.top().first < best_score) {
					best_cluster = i;
					best_score = c.pq.top().first;
				}
			}
		}

		int i = best_cluster;
		Cluster& c1a = clusters[i];
		int j = c1a.pq.top().second;
		Cluster& c2a = clusters[j];


		clusters.emplace_back(i, j, c1a.sz + c2a.sz);

		Cluster& c1 = clusters[i];
		Cluster& c2 = clusters[j];
		c1.pq.pop();
		c1.exists = 0;
		c2.exists = 0;




		for (int i = 0; i < clusters.size() - 1; i++) {
			if (! clusters[i].exists)
				continue;
			double score = 0;

			if (c1.d.count(i) && c2.d.count(i))
				score = (c1.sz * c1.d[i] + c2.sz * c2.d[i])/(c1.sz + c2.sz);
			else if (c1.d.count(i))
				score = c1.d[i];
			else if (c2.d.count(i))
				score = c2.d[i];
			else {

				continue;
			}
			clusters.back().pq.emplace(score, i);
			clusters.back().d[i] = score;
			clusters[i].d[clusters.size() - 1] = score;
		}

	}

	return clusters.back().newick(ts, clusters).str();

}

#endif
