#ifndef PROJECT_HPP
#define PROJECT_HPP

#include "DFN.hpp"
#include <vector>
#include <string>
#include <limits>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

bool readDFN(const string& filename, Fracture& frattura);
//bool printFracture(const Fracture& frattura);

vector<Polygon> createPolygons(const Fracture& frattura);

Plane calculatePlaneEquation(const Polygon& quad);

Vector3d calculateIntersectionDirection(const Plane& plane1, const Plane& plane2);

Vector3d calculateIntersectionPoint(const Plane& plane1, const Plane& plane2, const Vector3d& direction);
//bool printIntersection(const Plane& plane1, const Plane& plane2, int plane1Index, int plane2Index);

IntersectionLine calculateIntersectionLine(const Plane& plane1, const Plane& plane2);

string calculateLineEquation(const Point& p1, const Point& p2);

VerticesLine extractLinesFromPolygons(const vector<Polygon>& polygons);

Vector3d calculateLineIntersection(const Point& p1, const Vector3d& direction1, const Point& p2, const Vector3d& direction2);

bool isPointOnSegment(const Point& p1, const Point& p2, const Vector3d& intersection);

bool doSegmentsOverlap(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D, Vector3d& overlapStart, Vector3d& overlapEnd);

void calculateAndPrintIntersections(const vector<Polygon>& polygons, const IntersectionLine& intersectionLine, size_t i, size_t j, Traces& traces);

void saveTracesToFile(const string& filename, const Traces& traces);

bool isPointOnEdge(const Point& p1, const Point& p2, const Point& p);

double calculateDistance(const Point& p1, const Point& p2);

vector<TraceResult> checkTracePoints(const Traces& traces, const vector<Polygon>& polygons);

bool compareByLength(const TraceResult& a, const TraceResult& b);

void exportTraceResults(const vector<TraceResult>& results, const string& filename);

#endif // PROJECT_HPP
