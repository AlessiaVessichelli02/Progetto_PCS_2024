#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <iomanip>
#include <cmath>
#include <map>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct Fracture {

    unsigned int NumFractures = 0;
    vector<unsigned int> FractureID = {}; ///< ID delle fratture
    vector<unsigned int> NumVertices = {}; ///< Numero di vertici
    vector<MatrixXd> Vertices = {}; ///< Coordinate dei vertici vengono salvati in una matrice

};

struct Point {

    double x, y, z;

    Point() : x(0), y(0), z(0) {}

    Point(double x, double y, double z) : x(x), y(y), z(z) {}

};

struct Polygon {

    vector<Point> vertices; // Vettore contenente i vertici del poligono

};

struct Plane {

    double a, b, c, d; // Coefficienti dell'equazione del piano; ax + by + cz + d = 0
    vector<Point> vertices; // Per mantenere i vertici del piano

};

struct IntersectionLine {

    Vector3d direction; // Direzione della linea di intersezione
    Vector3d point; // Un punto sulla linea di intersezione

};

struct VerticesLine {

    vector<vector<string>> VerticesLines;

};

struct Traces {

    map<vector<int>,vector<vector<double>>> traces;
    map<int, vector<Point>> intersectionPoints;

};

// Define a structure to store trace data
struct TraceResult {

    int fractureId;
    int traceId;
    bool isNonPassante; // True for non-passant, False for passant
    double distance;

};

#endif // DFN_HPP
