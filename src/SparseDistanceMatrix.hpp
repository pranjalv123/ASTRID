#ifndef ASTRID_SPARSEDISTANCEMATRIX_ONCE__
#define ASTRID_SPARSEDISTANCEMATRIX_ONCE__

#include "phylokit/TaxonSet.hpp"
#include "phylokit/DistanceMatrix.hpp"

#include <unordered_map>
#include <utility>

class SparseDistanceMatrix {
    
private:
    struct hash_pair { 
        template <class T1, class T2> 
        size_t operator()(const std::pair<T1, T2>& p) const
        { 
            auto hash1 = std::hash<T1>{}(p.first); 
            auto hash2 = std::hash<T2>{}(p.second); 
            return hash1 ^ hash2; 
        } 
    }; 
    std::unordered_map<std::pair<Taxon, Taxon>, double, hash_pair> d;
    const TaxonSet& ts;
public:
    SparseDistanceMatrix(const TaxonSet& ts);
    SparseDistanceMatrix(const TaxonSet& ts, DistanceMatrix& distanceMatrix);
    SparseDistanceMatrix(const TaxonSet& ts, const std::string& newick);
    SparseDistanceMatrix();


    void set(Taxon t1, Taxon t2, double dist);
    double get(Taxon t1, Taxon t2) const;
    bool contains(Taxon t1, Taxon t2) const;

    std::unordered_map<std::pair<Taxon, Taxon>, double, hash_pair>::iterator begin() {return d.begin();}
    std::unordered_map<std::pair<Taxon, Taxon>, double, hash_pair>::iterator end() {return d.end();}
    std::unordered_map<std::pair<Taxon, Taxon>, double, hash_pair>::const_iterator begin() const noexcept {return d.begin();}
    std::unordered_map<std::pair<Taxon, Taxon>, double, hash_pair>::const_iterator end() const noexcept {return d.end();}
    std::unordered_map<std::pair<Taxon, Taxon>, double, hash_pair>::const_iterator cbegin() const noexcept {return d.cbegin();}
    std::unordered_map<std::pair<Taxon, Taxon>, double, hash_pair>::const_iterator cend() const noexcept {return d.cend();}

    DistanceMatrix toDistanceMatrix() const;
};

#endif