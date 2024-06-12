#include "src/Project.hpp"
#include "src/DFN.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;
//using namespace DFNLibrary;

int main() {

    Fracture frattura;

    if (!readDFN("FR3_data.txt", frattura)) {
        cerr << "Error" << endl;
        return 1;
    }

    vector<Polygon> polygons = createPolygons(frattura);
    VerticesLine verticesLine = extractLinesFromPolygons(polygons);

    cout << "Contenuto della struttura VerticesLine:" << endl;
    for (const auto& lines : verticesLine.VerticesLines) {
        for (const auto& line : lines) {
            cout << line << endl;
        }
        cout << "\n";
    }
    cout << "\n";

    vector<Plane> planes;
    for (const auto& poly : polygons) {
        planes.push_back(calculatePlaneEquation(poly));
    }

    for (size_t i = 0; i < planes.size(); ++i) {
        for (size_t j = i + 1; j < planes.size(); ++j) {
            IntersectionLine intersectionLine = calculateIntersectionLine(planes[i], planes[j]);
            cout << "Linea di intersezione tra il piano " << i << " e il piano " << j << ":" << endl;
            cout << "Punto sulla linea di intersezione: (" << intersectionLine.point.x() << ", " << intersectionLine.point.y() << ", " << intersectionLine.point.z() << ")" << endl;
            cout << "Vettore direzione: (" << intersectionLine.direction.x() << ", " << intersectionLine.direction.y() << ", " << intersectionLine.direction.z() << ")" << endl;
            cout << "Equazione della linea di intersezione: "
                 << "r(t) = (" << intersectionLine.point.x() << ", " << intersectionLine.point.y() << ", " << intersectionLine.point.z() << ") + t * ("
                 << intersectionLine.direction.x() << ", " << intersectionLine.direction.y() << ", " << intersectionLine.direction.z() << ")" << endl << endl;

            for (size_t k = 0; k < polygons.size(); ++k) {
                for (size_t l = 0; l < polygons[k].vertices.size(); ++l) {
                    Point p1 = polygons[k].vertices[l];
                    Point p2 = polygons[k].vertices[(l + 1) % polygons[k].vertices.size()];
                    Vector3d direction = Vector3d(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z).normalized();

                    try {
                        Vector3d intersection = calculateLineIntersection(p1, direction, {intersectionLine.point.x(), intersectionLine.point.y(), intersectionLine.point.z()}, intersectionLine.direction);
                        cout << "Punto di intersezione tra la retta del poligono " << k << " (vertici " << l << " e " << (l + 1) % polygons[k].vertices.size() << ") e la linea di intersezione tra i piani " << i << " e " << j << ": ("
                             << intersection.x() << ", " << intersection.y() << ", " << intersection.z() << ")" << endl;
                    } catch (const exception& e) {
                        cerr << "Errore nel calcolo dell'intersezione: " << e.what() << endl;
                    }
                }
            }
        }
    }

    Traces traces;
    // Stampa la struttura Traces
    cout << "Traces:" << endl;
    for (const auto& trace : traces.traces) {
        cout << "Fracture IDs: ";
        for (int fractureID : trace.first) {
            cout << fractureID << " ";
        }
        cout << endl;
        cout << "Intersection Points: ";
        for (const auto& point : trace.second) {
            cout << "(";
            for (int i = 0; i < point.size(); ++i) {
                cout << point[i];
                if (i < point.size() - 1) {
                    cout << ", ";
                }
            }
            cout << ")";
        }
        cout << endl;
    }

    saveTracesToFile("Traces_output.txt", traces);

    //checkTracePoints(traces, polygons);
    vector<TraceResult> results = checkTracePoints(traces, polygons);

    exportTraceResults(results, "TraceOutput.txt");

    return 0;
}


