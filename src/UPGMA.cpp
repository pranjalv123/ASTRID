#include "DistanceMethods.hpp"
#include "util/Logger.hpp"


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

		if (best_cluster == -1) { 
		  WARN << "Extremely high levels of missing data: species graph is disconnected.\nAccuracy will suffer.\n";
		  for (int i = 0; i < clusters.size(); i++) {
		    if (clusters[i].exists) {
		      best_cluster = clusters.size() - 1;
		      clusters.back().pq.emplace(9999, best_cluster);
		      break;
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
