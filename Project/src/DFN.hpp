#pragma once
#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace DFNLibrary{

struct DFN
{
    //Tracce
    unsigned int NumTraces = 0; ///< Numero di tracce
    unsigned int FractureID = 0; ///< ID delle fratture
    unsigned int NumVertices = 0; ///< Numero di vertici
    std::vector<Vector2i> Vertices = {}; ///< Vertici salvati in un vettore
};
}
