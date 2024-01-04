//
// Created by Pawulon on 31/12/2023.
//

#include "../Solver/Solver.hpp"


namespace mes::solver {

//    bool Solver::isOnEdge(int k) {
//        for (int i = 0; i <( width + 2 * nDampLayers); i++) {
//            if (i == 0 || i == (width + 2 * nDampLayers - 1)) {
//                for (int j = 0; j < height + 2 * nDampLayers; j++)
//                    if (k == j * width + 2 * nDampLayers + i) return true;
//
//            } else {
//                if (k == i || k == width * (height - 1) + i) return true;
//            }
//        }
//        return false;
//    }

    bool Solver::isOnEdge(int &i, int &j) {
        return i == 0 || j == 0 || i == (nVerWidth + 2 * nDampLayers - 1)|| j == (nVerHeight + 2 * nDampLayers - 1);
    }

   Solver& Solver::generateSimpleMesh() {
        for(int i = 0; i <  nVerWidth + 2 * nDampLayers; i++){
            for(int j = 0; j < nVerHeight + 2 * nDampLayers ; j++){
                points.emplace_back(-width * .5 - (nDampLayers * widthEleLen) + (widthEleLen * i),
                                    -height * .5 - (nDampLayers * heightEleLen) + (heightEleLen * j), isOnEdge(i,j));
            }
        }
//        for (int k = 0; k < (nVerWidth + 2 * nDampLayers) * (nVerHeight + 2 * nDampLayers); k++) {
//            int i = k % (nVerWidth + 2 * nDampLayers);
//            int j = k / (nVerWidth + 2 * nDampLayers);
//
//            points.emplace_back(-width / 2 - (nDampLayers * widthEleLen) + (widthEleLen * i),
//                                -height / 2 - (nDampLayers * heightEleLen) + (heightEleLen * j), isOnEdge(k));
//        }
        return initDampWalls();
    }

    Solver& Solver::initDampWalls() {
        if(nDampLayers>0){
            walls.insert(walls.begin(),   Wall::createFromDimensions(-width / 2 - (nDampLayers * widthEleLen), -height / 2 - (nDampLayers * heightEleLen),
                                          nDampLayers * widthEleLen, height + 2. * nDampLayers * heightEleLen,
                                           mes::DAMP));

            walls.insert(walls.begin(), Wall::createFromDimensions(-width / 2 - (nDampLayers * widthEleLen), -height / 2 - (nDampLayers * heightEleLen),
                                        width + 2. * nDampLayers * widthEleLen, nDampLayers * heightEleLen,
                                        mes::DAMP));

            walls.insert(walls.begin(), Wall::createFromDimensions(-width / 2 - (nDampLayers * widthEleLen), height / 2.,
                                         width + 2. * nDampLayers * widthEleLen, nDampLayers * heightEleLen,
                                         mes::DAMP));

            walls.insert(walls.begin(), Wall::createFromDimensions(width / 2., -height / 2 - nDampLayers * heightEleLen,
                                         nDampLayers * widthEleLen, height + 2. * nDampLayers * heightEleLen,
                                         mes::DAMP));
        }
        return *this;
    }

    Solver& Solver::setNumberOfDampLayers(int i) {
        this->nDampLayers = i;
        return *this;
    }

    Solver& Solver::divideIntoElements() {
        // with points rotation implemented inside
        fortunes::build(points, elementsIdx, walls, fType);

        // quad elements intersection

        nVertices = points.size();

        if (fType == mes::QUAD) {
            addQuadVertElements(points, elementsIdx);
            nVerticesQuad = static_cast<int>(points.size()) - nVertices;
            nVertices = nVertices + nVerticesQuad;
        }


        return *this;
    }

    Solver& Solver::buildStiffnessMatrix() {
        std::vector<Eigen::Triplet<std::complex<double>>> tripletListS{};

        for (auto &ele: Elements) {
            std::vector<std::vector<std::complex<double>>> &E = ele.E_;
            std::vector<int> &glob = ele.globalVectorIdx.indices;
            for (int k = 0; k < glob.size(); k++) {
                for (int l = 0; l < glob.size(); l++) {
                    int i = glob[k];
                    int j = glob[l];
                    if (!points[i].isBorder) {
                        tripletListS.emplace_back(i, j, E[k][l]);
                    }
                }
            }
        }

        for (int i = 0; i < points.size(); i++) {
            if (points[i].isBorder) {
                tripletListS.emplace_back(i, i, 1.);
            }
        }

        stiffnessMatrix = Eigen::SparseMatrix<std::complex<double>>(nVertices, nVertices);
        stiffnessMatrix.setFromTriplets(tripletListS.begin(), tripletListS.end());

        return *this;
    }

    Solver& Solver::buildLoadVector() {
        loadVector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>(nVertices);
        int minIdx = 0;
        double minDist = sqrt(pow(points[minIdx].x - (sourcePoint.x), 2) + pow(points[minIdx].y - (sourcePoint.y), 2));
        loadVector[0] = 0;
        for (int i = 0; i < nVertices - nVerticesQuad; i++) {
            if(!points[i].isBorder){
                double temp = sqrt(pow(points[i].x - (sourcePoint.x), 2) + pow(points[i].y - (sourcePoint.y), 2));
                if (temp < minDist) {
                    minDist = temp;
                    minIdx = i;
                }
            }
            loadVector[i] = 0;
        }
        loadVector[minIdx] = 1.;
        return *this;
    }

    Solver& Solver::buildSolver() {

        solverLU.compute(stiffnessMatrix);
        return *this;
    }

    Solver& Solver::solve() {
        buildLoadVector();
        solutions = solverLU.solve(loadVector);

        maxSolutionValue = 0.;
        for (auto &v: solutions) {
            if (abs(v) > maxSolutionValue) {
                maxSolutionValue = abs(v);
            }
        }

        return *this;
    }

    Solver& Solver::solve(double sx, double sy) {
        setSourcePoint(sx, sy);
        return solve();
    }

    Solver& Solver::draw() {
        auto normalSolution = normalizeSolution(); //placeholder
        plot::CreateOpenGlWindow(points, Elements, normalSolution, width, height, canvasWidth, canvasHeight);
        return *this;
    }

    std::vector<double> Solver::normalizeSolution() {
        std::vector<double> out;

        for (auto &v: solutions) {
            out.push_back(abs(v) / maxSolutionValue);
        }

        return out;
    }


    Solver::Solver(double width, double height, int nVerWidth, int nVerHeight) :
                                                                            width(width), height(height),
                                                                            nVerWidth(nVerWidth),
                                                                            nVerHeight(nVerHeight) {
        widthEleLen = width / (nVerWidth - 1);
        heightEleLen = height / (nVerHeight - 1);
        canvasWidth = 900;
        canvasHeight = static_cast<unsigned int>(canvasWidth * (height/width));

        mes::TriangleElement::k_ = 2. * M_PI * frequency * pow(10, 9) / (2.99792458 * pow(10, 8));
    }

    Solver& Solver::buildElements() {
        int m = 0;
        for (auto &idx: elementsIdx) {
            Elements.emplace_back(m++, points, idx);
        }
        return *this;
    }

    Solver &Solver::addWall(Wall w) {
        walls.push_back(w);
        return *this;;
    }

    Solver &Solver::setSourcePoint(double x, double y) {
        this->sourcePoint = Vertex2D(x,y);
        return *this;
    }


    Solver &Solver::setFrequency(double f) {
        this->frequency = f;
        mes::TriangleElement::k_ = 2. * M_PI * frequency * pow(10, 9) / (2.99792458 * pow(10, 8));
        return *this;
    }

    Solver &Solver::setBaseFunctionType(mes::baseFuncType fType) {
        this->fType = fType;
        return *this;
    }

    Solver &Solver::setImageSize(unsigned int x, unsigned int y) {
        canvasWidth = x;
        canvasHeight = y;
        return *this;
    }

    Solver &Solver::buildStructure(bool compTime) {
        if(compTime) {
            MEASURE_TIME(generateSimpleMesh());
            MEASURE_TIME(divideIntoElements());
            MEASURE_TIME(buildElements());
            return *this;
        }
        return generateSimpleMesh().divideIntoElements().buildElements();
    }

    Solver &Solver::computeSolver(bool compTime) {
        if(compTime) {
            MEASURE_TIME(buildStiffnessMatrix());
            MEASURE_TIME(buildSolver());
            return *this;
        }
        return buildStiffnessMatrix().buildSolver();
    }

    double Solver::getMaxValue() const {
        return maxSolutionValue;
    }


}


