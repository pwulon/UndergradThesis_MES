//
// Created by Pawulon on 28/12/2023.
//

#include "Section.hpp"


struct Section{
    int i1_, i2_, i3_;

    Section(int _i1, int _i2, int _i3 = 0): i3_{_i3}{
        if(_i1 < _i2 ){
            i1_ = _i1;
            i2_ = _i2;
        }else{
            i1_ = _i2;
            i2_ = _i1;
        }
    }
};

struct SecComp{
    constexpr bool operator()(const Section &p1, const Section &p2) const{
        return (p1.i1_ == p2.i1_  && p1.i2_  == p2.i2_);
    }
};

struct SecHash {
    std::size_t operator () (const Section& p) const {
        auto h1 = std::hash<int>{}(p.i1_);
        auto h2 = std::hash<int>{}(p.i2_);
//        auto h3 = std::hash<T1>{}(p[2]);

        // A simple combination hash function
        return h1 ^ h2;
    }
};
void addQuadVertElements(std::vector<Vertex2D> &ver, std::vector<std::vector<int>> elementsIdx){
//TODO bug dla isBorder, nie wiem jak zdecydowac czy nowy punkt jest na brzegu dla warunkow brzegowych direchleta;

    std::unordered_set<Section, SecHash, SecComp> test_set;
    for(auto &idxs:elementsIdx){
        for(int i=0;i<3;i++){
            auto iter = test_set.find({idxs[i],idxs[(i+1)%3]});
            if ( iter != test_set.end()){
                idxs.push_back(iter->i3_);
            }else{
                ver.emplace_back((ver[idxs[i]].x + ver[idxs[(i+1)%3]].x)/2, (ver[idxs[i]].y + ver[idxs[(i+1)%3]].y)/2);
//                ver.emplace_back(v, ver[idxs[i]].isBorder && ver[idxs[(i+1)%3]].isBorder);
                int i3 = static_cast<int>(ver.size() - 1);
                test_set.insert({idxs[i],idxs[(i+1)%3],i3});
                idxs.push_back(i3);
            }
        }
    }
}