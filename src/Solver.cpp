//
// Created by Pawulon on 31/12/2023.
//

#include "../Solver/Solver.hpp"


namespace mes::solver {

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
        return initDampWalls();
    }

    Solver& Solver::initDampWalls() {
        if(nDampLayers>0){
            walls.insert(walls.begin(),  Wall::createFromDimensions(-width / 2 - (nDampLayers * widthEleLen), -height / 2 - (nDampLayers * heightEleLen),
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


    Solver &Solver::setDampLayersDepth(double d) {
        this->nDampLayers = static_cast<int>(d/dx);
        return *this;
    }

    Solver& Solver::divideIntoElements() {
        // with points rotation implemented inside

        fortunes::build(points, elementsIdx, walls, fType);

        // quad elements intersection

        nVertices = nVerticesLin = static_cast<int>(points.size());

        if (fType == mes::QUAD) {
            addQuadVertElements(points, elementsIdx);
            nVerticesLin = nVertices;
            nVertices = static_cast<int>(points.size());
        }


        for(auto &ei:elementsIdx){
            if(ei.et != mes::AIR){
                for(auto &i: ei.indices) {
                    points[i].setElementType(ei.et);
                }
            }
        }

        return *this;
    }

    Solver& Solver::buildElements() {
        for (auto &idx: elementsIdx) {
            Elements.emplace_back(points, idx);
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
                tripletListS.emplace_back(i, i, 1.f);
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
        loadVector[0] = 0.f;
        for (int i = 0; i < nVertices; i++) {
            if(i < nVerticesLin && !points[i].isBorder){
                double temp = sqrt(pow(points[i].x - (sourcePoint.x), 2) + pow(points[i].y - (sourcePoint.y), 2));
                if (temp < minDist) {
                    minDist = temp;
                    minIdx = i;
                }
            }
            loadVector[i] = 0.f;
        }
        loadVector[minIdx] = 1.f;
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

    std::string Solver::draw() {
        auto normalSolution = normalizeSolutionAbs(); //placeholder
        auto filename = plot::CreateOpenGlWindow(points, Elements, normalSolution, width, height, canvasWidth, canvasHeight, this->fType, true);
        return filename;
    }

    Solver& Solver::drawImag() {
        auto normalSolution = normalizeSolutionImg(); //placeholder
        plot::CreateOpenGlWindow(points, Elements, normalSolution, width, height, canvasWidth, canvasHeight,this->fType, false);
        return *this;
    }

    Solver& Solver::drawReal() {
        auto normalSolution = normalizeSolutionReal(); //placeholder
        plot::CreateOpenGlWindow(points, Elements, normalSolution, width, height, canvasWidth, canvasHeight, this->fType, false);
        return *this;
    }

    std::vector<double> Solver::normalizeSolutionAbs() {
        std::vector<double> out;

        for (auto &v: solutions) {
            out.push_back(abs(v) / maxSolutionValue);
        }

        return out;
    }

    std::vector<double> Solver::normalizeSolutionImg() {
        std::vector<double> out(nVertices);

        double maxSolutionValueReal = 0.;
        for (auto &v: solutions) {
            if (v.imag() > maxSolutionValueReal) {
                maxSolutionValueReal = v.imag();
            }
        }

        for (int i=0;i<nVertices;i++) {
            out[i] = ((solutions[i].imag() +  maxSolutionValueReal) / (2.*maxSolutionValueReal));
        }

        return out;
    }

    std::vector<double> Solver::normalizeSolutionReal() {
        std::vector<double> out;

        double maxSolutionValueReal = 0.;
        for (auto &v: solutions) {
            if (v.real() > maxSolutionValueReal) {
                maxSolutionValueReal = v.real();
            }
        }

        for (auto &v: solutions) {
            out.push_back((v.real() +  maxSolutionValueReal) / (2.*maxSolutionValueReal));
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

        dx = std::min(widthEleLen, heightEleLen);
        setDampLayersDepth(.25);

        mes::TriangleElement::k_ = 2. * M_PI * frequency * pow(10, 9) / (2.99792458 * pow(10, 8));
    }

    Solver::Solver(double width, double height, double _dx) :dx{_dx},
                                             width(width), height(height)
    {
        nVerWidth = static_cast<int>(width/dx) + 1;
        nVerHeight = static_cast<int>(height/dx) + 1;

        widthEleLen = heightEleLen = dx;
        canvasWidth = 900;
        canvasHeight = static_cast<unsigned int>(canvasWidth * (height/width));
        setDampLayersDepth(.25);

        mes::TriangleElement::k_ = 2. * M_PI * frequency * pow(10, 9) / (2.99792458 * pow(10, 8));
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


