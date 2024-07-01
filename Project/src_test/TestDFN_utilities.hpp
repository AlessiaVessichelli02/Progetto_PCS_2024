#ifndef __UCDUtilities_H
#define __UCDUtilities_H

#include "Eigen/Eigen"
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;
using namespace Eigen;

struct Point {

    double x, y, z;

    Point() : x(0), y(0), z(0) {}

    Point(double x, double y, double z) : x(x), y(y), z(z) {}

};

struct IntersectionLine {

    Vector3d direction; // Direzione della linea di intersezione
    Vector3d point; // Un punto sulla linea di intersezione

};

struct Plane {

    double a, b, c, d; // Coefficienti dell'equazione del piano; ax + by + cz + d = 0
    vector<Point> vertices; // Per mantenere i vertici del piano

};

struct Polygon {

    vector<Point> vertices; // Vettore contenente i vertici del poligono

};

// Funzione per calcolare la distanza euclidea tra due punti
double calculateDistance(const Point& p1, const Point& p2);

bool isPointOnSegment(const Point& p1, const Point& p2, const Vector3d& intersection);

Vector3d calculateIntersectionDirection(const Plane& plane1, const Plane& plane2);

Vector3d calculateIntersectionPoint(const Plane& plane1, const Plane& plane2, const Vector3d& direction);

Vector3d calculateIntersectionBetweenLines(const Point& p1, const Vector3d& direction1, const Point& p2, const Vector3d& direction2);

bool doSegmentsOverlap(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D, Vector3d& overlapStart, Vector3d& overlapEnd);

double calculateSphereRadius(const Polygon& poly);

bool doPolygonsIntersect(const Polygon& poly1, const Polygon& poly2);


#endif

