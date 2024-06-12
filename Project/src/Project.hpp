#ifndef PROJECT_HPP
#define PROJECT_HPP

#include "DFN.hpp"
#include <vector>
#include <string>

using namespace std;

//namespace DFNLibrary {

struct DFN {

    vector<Fracture> fractures;

    void readDFN(const string& filename);
    void print() const;
};

#endif // PROJECT_HPP

//} // namespace DFN



/*double calculateLength(const Point3D& p1, const Point3D& p2);

vector<Fracture> readDFNFromFile(const string& filename);

void writeTracesToFile(const string& filename, const vector<Trace>& traces);

void writeFractureTracesToFile(const string& filename, const vector<vector<Trace>>& fractureTraces);

bool doLinesIntersect(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4, Point3D& intersection);

bool isPointOnLineSegment(const Point3D& p, const Point3D& p1, const Point3D& p2);
*/
