//
// Created by Pawulon on 31/12/2023.
//

#include "Solver.hpp"


bool fem::solve::Solver::isOnEdge(int k){
    for(int i = 0; i<width; i++){
        if(i==0 || i==width - 1){
            for(int j=0;j<height;j++)
                if(k == j*width + i) return true;

        }else
        {
            if(k == i || k == width*(height - 1) + i) return true;
        }
    }
    return false;
}

void fem::solve::Solver::generateSimpleMesh(){


    for(int k =0; k<(nVerWidth)*(nVerHeight);k++){
        int i = k%(nVerWidth);
        int j = k/(nVerWidth);

//        points.emplace_back(- width/2 + (widthEleLen * i), - height/2 + (heightEleLen * j ), isBorder(k, nVerWidth, nVerHeight));
    }
}