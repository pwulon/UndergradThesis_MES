//
// Created by Pawulon on 31/12/2023.
//

#include "../Solver/Solver.hpp"


namespace mes::solver {

    bool Solver::isOnEdge(int &i, int &j) {
        return i == 0 || j == 0 || i == (nVerWidth + 2 * nDampLayers - 1)|| j == (nVerHeight + 2 * nDampLayers - 1);
    }

    Solver::Solver(mes::baseFuncType _bft, double _f): frequency{_f}, fType{_bft} {

        TriangleElement::staticMembersInit(this->fType);


        canvasWidth = 1800;
        canvasHeight = 1800;
        mes::setRefIdx(frequency);
        mes::TriangleElement::k_ = 2. * M_PI * frequency * 10. / 2.99792458 ;
    }

   Solver& Solver::generateSimpleMesh() {

//       int numPoints=50;
//       std::srand(static_cast<unsigned int>(std::time(nullptr)));
//
//       // Generate and print n random points in the 1 by 1 square
//       for (int i = 0; i < numPoints; ++i) {
//           double x = static_cast<double>(std::rand()) / RAND_MAX;  // Random x-coordinate between 0 and 1
//           double y = static_cast<double>(std::rand()) / RAND_MAX;  // Random y-coordinate between 0 and 1
//
////           std::cout << "(" << x << ", " << y << ")" << std::endl;
//           points.emplace_back(x, y);
//       }



//       points.emplace_back(0., 0.);
//       for (int i = 1; i < 5; ++i) {
//           double radius = i * dx;
//           for (int j = 0; j < numPoints; ++j) {
//               double theta = 2.0 * M_PI * j / numPoints;
//               double x = radius * std::cos(theta);
//               double y = radius * std::sin(theta);
//
//               points.emplace_back(x, y, i==49);
////               std::cout << "(" << x << ", " << y << ")" << std::endl;
//           }
//           numPoints+=10;
//       }

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
                                           mes::CONCRETE));

            walls.insert(walls.begin(), Wall::createFromDimensions(-width / 2 - (nDampLayers * widthEleLen), -height / 2 - (nDampLayers * heightEleLen),
                                        width + 2. * nDampLayers * widthEleLen, nDampLayers * heightEleLen,
                                        mes::CONCRETE));

            walls.insert(walls.begin(), Wall::createFromDimensions(-width / 2 - (nDampLayers * widthEleLen), height / 2.,
                                         width + 2. * nDampLayers * widthEleLen, nDampLayers * heightEleLen,
                                         mes::CONCRETE));

            walls.insert(walls.begin(), Wall::createFromDimensions(width / 2., -height / 2 - nDampLayers * heightEleLen,
                                         nDampLayers * widthEleLen, height + 2. * nDampLayers * heightEleLen,
                                         mes::CONCRETE));
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

//        std::cout<<points.size()<<std::endl;
        if (fType == mes::QUAD) {
            addQuadVertElements(points, elementsIdx);
            nVerticesLin = nVertices;
            nVertices = static_cast<int>(points.size());
        }
//        std::cout<<points.size()<<std::endl;
//        std::cout<<elementsIdx.size()<<std::endl;
//        for(auto &ei:elementsIdx){
//            if(ei.et != mes::AIR){
//                for(auto &i: ei.indices) {
//                    points[i].setElementType(ei.et);
//                }
//            }
//        }

//        std::ofstream filep, filee;
//        filep.open("point.txt");
//        filee.open("elements.txt");
//
//        for(auto&p:points){
//            filep<<p.x<<" "<<p.y<<std::endl;
//        }
//        for(auto & e:elementsIdx){
//
//            for(auto& i: e.indices){
//                filee<<i<<" ";
//            }
//            filee<<"\n";
//        }
//        filep.close();
//        filee.close();

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

        minIdx = 0;
        double minDist = sqrt(pow(points[minIdx].x - (sourcePoint.x), 2) + pow(points[minIdx].y - (sourcePoint.y), 2));
//        loadVector[0] = 0.f;
        for (int i = 0; i < nVertices; i++) {
            if(i < nVerticesLin && !points[i].isBorder){
                double temp = sqrt(pow(points[i].x - (sourcePoint.x), 2) + pow(points[i].y - (sourcePoint.y), 2));
                if (temp < minDist) {
                    minDist = temp;
                    minIdx = i;
                }
            }
//            loadVector[i] = 0.f;
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
        return *this;
    }

    Solver& Solver::solve(double sx, double sy, bool add) {

        setSourcePoint(sx, sy).solve();

        if(add ){
            for(int i=0;i<nVertices;i++){
                sum_solutions[i] = abs(solutions[i]);
            }
        }else{
            loadVector[minIdx] = 0.;
            for(int i=0;i<nVertices;i++){
                sum_solutions[i] = abs(solutions[i]);
            }
        }


        return *this;
    }

    std::string Solver::draw() {
        auto normalSolution = normalizeSolutionAbs(); //placeholder
        auto filename = plot::CreateOpenGlWindow(points, Elements, normalSolution, width, height, canvasWidth, canvasHeight, this->fType, true);
        return filename;
    }

    std::string Solver::drawImag() {
        auto normalSolution = normalizeSolutionImg(); //placeholder
        auto filename = plot::CreateOpenGlWindow(points, Elements, normalSolution, width, height, canvasWidth, canvasHeight, this->fType, false);
        return filename;
    }

    std::string Solver::drawReal() {
        auto normalSolution = normalizeSolutionReal(); //placeholder
        auto filename = plot::CreateOpenGlWindow(points, Elements, normalSolution, width, height, canvasWidth, canvasHeight, this->fType, false);
        return filename;
    }

    std::vector<double> Solver::normalizeSolutionAbs() {
        std::vector<double> out;

        maxSolutionValue = 0.;
        for (auto &v: sum_solutions) {
            if (abs(v) > maxSolutionValue) {
                maxSolutionValue = abs(v);
            }
        }

        for (auto &v: sum_solutions) {
            out.push_back(abs(v) / maxSolutionValue);
        }

        return out;
    }

    std::vector<double> Solver::normalizeSolutionImg() {
        std::vector<double> out;

        double maxSolutionValueReal = 0.;
        for (auto &v: solutions) {
            if (abs(v.imag()) > maxSolutionValueReal) {
                maxSolutionValueReal = abs(v.imag());
            }
        }

        for (auto &v: solutions) {
            out.push_back((v.imag() +  maxSolutionValueReal) / (2.*maxSolutionValueReal));
        }

        maxSolutionValue = maxSolutionValueReal;

        return out;
    }

    std::vector<double> Solver::normalizeSolutionReal() {
        std::vector<double> out;

        double maxSolutionValueReal = 0.;
        for (auto &v: solutions) {
            if (abs(v.real()) > maxSolutionValueReal) {
                maxSolutionValueReal = abs(v.real());
            }
        }

        for (auto &v: solutions) {
            out.push_back((v.real() +  maxSolutionValueReal) / (2.*maxSolutionValueReal));
        }
        maxSolutionValue = maxSolutionValueReal;

        return out;
    }


    Solver& Solver::setMeshSize(double _width, double _height, int _nVerWidth, int _nVerHeight) {
        width = _width;
        height = _height;
        nVerWidth = _nVerWidth;
        nVerHeight = _nVerHeight;

        dx = std::min(widthEleLen, heightEleLen);

        widthEleLen = width / (nVerWidth - 1);
        heightEleLen = height / (nVerHeight - 1);

        setDampLayersDepth(.25);

        return *this;
    }

    Solver& Solver::setMeshSize(double _width, double _height, double _dx)
    {
        width = _width;
        height = _height;
        dx = _dx;
        nVerWidth = static_cast<int>(width/dx) + 1;
        nVerHeight = static_cast<int>(height/dx) + 1;

        widthEleLen = heightEleLen = dx;
        setDampLayersDepth(.5);
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
        loadVector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>(nVertices);
        sum_solutions = std::vector<double>(nVertices, 0.);
        for (int i = 0; i < nVertices; i++) {
            loadVector[i] = 0.f;
        }


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


