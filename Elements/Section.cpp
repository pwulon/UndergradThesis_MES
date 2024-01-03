//
// Created by Pawulon on 28/12/2023.
//

#include "Section.hpp"


struct Section{
    int i1_, i2_, i3_;

    Section(int &_i1, int &_i2, int _i3 = 0): i3_{_i3}{
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
void addQuadVertElements(std::vector<Vertex2D> &ver, std::vector<fem::ElementIndices> &elementsIdx){
//TODO bug dla isBorder, nie wiem jak zdecydowac czy nowy punkt jest na brzegu dla warunkow brzegowych direchleta;

    std::unordered_set<Section, SecHash, SecComp> sectionSet;
    for(auto &eli:elementsIdx){
        auto &idxs = eli.indices;
        for(int i=0;i<3;i++){
            const Vertex2D& vertex1 = ver[idxs[i]];
            const Vertex2D& vertex2 = ver[idxs[(i + 1) % 3]];
            auto iter = sectionSet.find({idxs[i], idxs[(i + 1) % 3]});
            if (iter != sectionSet.end()){
                idxs.push_back(iter->i3_);
                sectionSet.erase(iter);
            }else{
                ver.emplace_back((vertex1.x + vertex2.x)/2, (vertex1.y + vertex2.y)/2,
                                 vertex1.isBorder && vertex2.isBorder);

                int i3 = static_cast<int>(ver.size() - 1);
                sectionSet.insert({idxs[i], idxs[(i + 1) % 3], i3});
                idxs.push_back(i3);
            }
        }
        //magic rearanging
        idxs = {idxs[5], idxs[0], idxs[3], idxs[1], idxs[4], idxs[2]};
    }
}