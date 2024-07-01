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

double calculateDistance(const Point& p1, const Point& p2);

double calculateCircumferenceRadius(const Polygon& poly);

bool doPolygonsIntersect(const Polygon& poly1, const Polygon& poly2);

Plane calculatePlaneEquation(const Polygon& quad);

Vector3d calculateIntersectionDirection(const Plane& plane1, const Plane& plane2);

Vector3d calculateIntersectionPoint(const Plane& plane1, const Plane& plane2, const Vector3d& direction);

IntersectionLine calculateIntersectionLine(const Plane& plane1, const Plane& plane2);

Vector3d calculateIntersectionBetweenLines(const Point& p1, const Vector3d& direction1, const Point& p2, const Vector3d& direction2);

bool isPointOnSegment(const Point& p1, const Point& p2, const Vector3d& intersection);

bool doSegmentsOverlap(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D, Vector3d& overlapStart, Vector3d& overlapEnd);

void calculateAndPrintIntersections(const vector<Polygon>& polygons, const IntersectionLine& intersectionLine, size_t i, size_t j, Traces& traces);

void saveTracesToFile(const string& filename, const Traces& traces);

bool compareByLength(const TraceResult& a, const TraceResult& b);

void checkTracePoints(const Traces& traces, const vector<Polygon>& polygons, TraceResult& traceResult);
//void printTraceResult(const TraceResult& traceResult);

void exportTraceResult(const string& filename, const TraceResult& traceResult);

#endif // PROJECT_HPP
