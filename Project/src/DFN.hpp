#ifndef DFN_HPP
#define DFN_HPP

#include <vector>
#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

struct Vertex {
    double x, y, z;
    Vertex(double x, double y, double z) : x(x), y(y), z(z) {}
};

struct Fracture {
    int id;
    vector<Vertex> vertices;

    Fracture(int id) : id(id) {}

    void addVertex(const Vertex& vertex) {
        vertices.push_back(vertex);
    }

    void print() const {
        cout << "Fracture ID: " << id << "\n";
        for (const auto& vertex : vertices) {
            std::cout << "Vertex: (" << setprecision(16) << vertex.x << ", " << setprecision(16) << vertex.y << ", " << setprecision(16) << vertex.z << ")\n";
        }
    }
};

#endif // DFN_HPP


/*
#include <iostream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace DFNLibrary{

struct DFN
{
    //Tracce
    unsigned int NumFractures = 0;
    unsigned int NumTraces = 0; ///< Numero di tracce
    unsigned int FractureID = 0; ///< ID delle fratture
    unsigned int NumVertices = 0; ///< Numero di vertici
    vector<Vector3d> CoordinatesVertices = {}; ///< Coordinate dei vertici vengono salvati in un vettore

};

}

#endif
*/

/*
#pragma once

#include <vector>
#include <string>

using namespace std;

namespace DFN {

struct Point3D {
    double x, y, z;
};

struct Fracture {
    int id;
    vector<Point3D> vertices;
};



struct Trace {
    int id;
    int fractureId1;
    int fractureId2;
    Point3D point1;
    Point3D point2;
    double length;
    bool tips;
};

struct MyDFN {
    vector<Fracture> fractures;
    vector<Trace> traces;
    vector<vector<Trace>> fractureTraces;

    MyDFN(const string& filename);
    void calculateTraces();
    void findTracesBetweenFractures(const Fracture& f1, const Fracture& f2);
    void classifyTrace(const Trace& trace, vector<Trace>& passTraces, vector<Trace>& nonPassTraces);
    void classifyAndSortTraces();
    void writeResults(const string& traceFilename, const string& fractureTraceFilename);
};

} // namespace DFN
*/
