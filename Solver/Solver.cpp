//
// Created by Pawulon on 31/12/2023.
//

#include "Solver.hpp"


namespace fem::solve {

    bool Solver::isOnEdge(int k) {
        for (int i = 0; i < width + 2 * nDampLayers; i++) {
            if (i == 0 || i == width + 2 * nDampLayers - 1) {
                for (int j = 0; j < height + 2 * nDampLayers; j++)
                    if (k == j * width + 2 * nDampLayers + i) return true;

            } else {
                if (k == i || k == width * (height - 1) + i) return true;
            }
        }
        return false;
    }

   Solver& Solver::generateSimpleMesh() {
        for (int k = 0; k < (nVerWidth + 2 * nDampLayers) * (nVerHeight + 2 * nDampLayers); k++) {
            int i = k % (nVerWidth + 2 * nDampLayers);
            int j = k / (nVerWidth + 2 * nDampLayers);

            points.emplace_back(-width / 2 - (nDampLayers * widthEleLen) + (widthEleLen * i),
                                -height / 2 - (nDampLayers * heightEleLen) + (heightEleLen * j), isOnEdge(k));
        }
        return initDampWalls();
    }

    Solver& Solver::initDampWalls() {
        walls.insert(walls.begin(),   {-width / 2 - (nDampLayers * widthEleLen), -height / 2 - (nDampLayers * heightEleLen),
                                      nDampLayers * widthEleLen, height + 2. * nDampLayers * heightEleLen,
                                      fem::DAMP});

        walls.insert(walls.begin(),{-width / 2 - (nDampLayers * widthEleLen), -height / 2 - (nDampLayers * heightEleLen),
                                    width + 2. * nDampLayers * widthEleLen, nDampLayers * heightEleLen,
                                    fem::DAMP});

        walls.insert(walls.begin(), {-width / 2 - (nDampLayers * widthEleLen), height / 2.,
                                     width + 2. * nDampLayers * widthEleLen, nDampLayers * heightEleLen,
                                     fem::DAMP});

        walls.insert(walls.begin(), {width / 2., -height / 2 - nDampLayers * heightEleLen,
                                     nDampLayers * widthEleLen, height + 2. * nDampLayers * heightEleLen,
                                     fem::DAMP});
        return *this;
    }

    Solver& Solver::setNumberOfDampLayers(int i) {
        this->nDampLayers = i;
        return *this;
    }

    Solver& Solver::divideIntoElements() {
        // with points rotation implemented inside
        build(points, elementsIdx, walls, fType);

        // quad elements intersection
        if (fType == fem::QUAD) {
            addQuadVertElements(points, elementsIdx);
        }

        nVertices = points.size();
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
        for (int i = 1; i < nVertices; i++) {
            double temp = sqrt(pow(points[i].x - (sourcePoint.x), 2) + pow(points[i].y - (sourcePoint.y), 2));
            if (temp < minDist) {
                minDist = temp;
                minIdx = i;
            }
            loadVector[i] = 0;
        }
        loadVector[minIdx] = 1.;
        return *this;
    }

    Solver& Solver::buildLoadVector(double sx, double sy) {
        return setSourcePoint(sx,sy).buildLoadVector();
    }

    Solver &Solver::buildLoadVector(const Vertex2D &p) {
        return setSourcePoint(p).buildLoadVector();
    }

    Solver& Solver::buildSolver() {
        solverLU.compute(stiffnessMatrix);
        return *this;
    }

    Solver& Solver::solve() {
        solutions = solverLU.solve(loadVector);
        return *this;
    }

    std::string Solver::draw() {
        auto normalSolution = normalizeSolution(); //placeholder
        auto encodedImg = CreateOpenGlWindow(points, Elements, normalSolution, width, height, canvasWidth, canvasHeight);
        return encodedImg;
    }

    std::vector<double> fem::solve::Solver::normalizeSolution() {
        std::vector<double> out;

        double maxElement = 0;
        for (auto &v: solutions) {
            if (abs(v) > maxElement) {
                maxElement = abs(v);
            }
        }
        for (auto &v: solutions) {
            out.push_back(abs(v) / maxElement);
        }

        std::cout << "abs result range: [" << 0 << ", " << maxElement << "]\n";
        return out;
    }

    Solver& Solver::doAll() {
        generateSimpleMesh()
            .initDampWalls()
            .divideIntoElements()
            .buildElements()
            .buildStiffnessMatrix()
            .buildSolver()
            .buildLoadVector()
            .solve()
            .draw();

        return *this;
    }

    Solver::Solver(double width, double height, int nVerWidth, int nVerHeight) :
                                                                            width(width), height(height),
                                                                            nVerWidth(nVerWidth),
                                                                            nVerHeight(nVerHeight) {
        widthEleLen = width / (nVerWidth - 1);
        heightEleLen = height / (nVerHeight - 1);
        canvasWidth = 900;
        canvasHeight = static_cast<unsigned int>(canvasWidth * (height/width));

        fem::TriangleElement::k_ = 2. * M_PI * frequency * pow(10,9)/(2.99792458 * pow(10,8));
    }

    Solver& Solver::buildElements() {
        int m = 0;
        auto start = std::chrono::high_resolution_clock::now();
        for (auto &idx: elementsIdx) {
            Elements.emplace_back(m++, points, idx);
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
        std::cout << "Time taken by function: " << duration.count()/Elements.size() << std::endl;

        return *this;
    }

    Solver &Solver::addWall(double leftDownX, double leftDownY, double w, double h, fem::elementType elt) {
        walls.emplace_back(leftDownX, leftDownY, w, h, elt);
        return *this;
    }

    Solver &Solver::addWall(Wall &w) {
        walls.push_back(w);
        return *this;;
    }

    Solver &Solver::setSourcePoint(double x, double y) {
        this->sourcePoint = Vertex2D(x,y);
        return *this;
    }

    Solver &Solver::setSourcePoint(const Vertex2D &p) {
        this->sourcePoint = Vertex2D(p);
        return *this;
    }

    Solver &Solver::setFrequency(double f) {
        this->frequency = f;
        fem::TriangleElement::k_ = 2. * M_PI * frequency * pow(10,9)/(2.99792458 * pow(10,8));
        return *this;
    }

    Solver &Solver::setBaseFunctionType(fem::baseFuncType fType) {
        this->fType = fType;
        return *this;
    }

    Solver &Solver::setImageSize(unsigned int x, unsigned int y) {
        canvasWidth = x;
        canvasHeight = y;
        return *this;
    }


}


