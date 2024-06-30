#include "src/Project.hpp"
#include "src/DFN.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace Eigen;

int main() {

    Fracture frattura;

    if (!readDFN("FR10_data.txt", frattura)) {
        cerr << "Error" << endl;
        return 1;
    }

    vector<Polygon> polygons = createPolygons(frattura);

    vector<Plane> planes;
    for (size_t i = 0; i < polygons.size(); ++i) {
        const auto& poly = polygons[i];
        cout << "Polygon " << i << " vertices:\n";
        for (size_t j = 0; j < poly.vertices.size(); ++j) {
            cout << "P" << j + 1 << ": (" << poly.vertices[j].x << ", " << poly.vertices[j].y << ", " << poly.vertices[j].z << ")\n";
        }
        Plane plane = calculatePlaneEquation(poly);
        planes.push_back(plane);
        cout << "Equation of the plane: " << plane.a << "x + " << plane.b << "y + " << plane.c << "z + " << plane.d << " = 0\n\n";
    }

    Traces traces; // Inizializzazione della struttura Traces

    // Itera su tutti i possibili coppie di poligoni
    for (size_t i = 0; i < polygons.size(); ++i) {
        for (size_t j = i + 1; j < polygons.size(); ++j) {

            // Se doPolygonsIntersect ritorna True, allora si può calcolare le intersezioni
            if (doPolygonsIntersect(polygons[i], polygons[j])) {

                // Calcola la linea di intersezione tra i piani dei poligoni i e j
                Plane plane1 = calculatePlaneEquation(polygons[i]);
                Plane plane2 = calculatePlaneEquation(polygons[j]);
                IntersectionLine intersectionLine = calculateIntersectionLine(plane1, plane2);

                // Stampa le informazioni sulla linea di intersezione
                cout << "Linea di intersezione tra il poligono " << i << " e il poligono " << j << ":" << endl;
                cout << "Punto sulla linea di intersezione: (" << intersectionLine.point.x() << ", " << intersectionLine.point.y() << ", " << intersectionLine.point.z() << ")" << endl;
                cout << "Vettore direzione: (" << intersectionLine.direction.x() << ", " << intersectionLine.direction.y() << ", " << intersectionLine.direction.z() << ")" << endl;
                cout << "Equazione della linea di intersezione: "
                     << "r(t) = (" << intersectionLine.point.x() << ", " << intersectionLine.point.y() << ", " << intersectionLine.point.z() << ") + t * (" << intersectionLine.direction.x() << ", " << intersectionLine.direction.y() << ", " << intersectionLine.direction.z() << ")" << endl;

                // Calcola e stampa le intersezioni effettive
                calculateAndPrintIntersections(polygons, intersectionLine, i, j, traces);
            }
        }
    }


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

    TraceResult traceResult;

    checkTracePoints(traces, polygons, traceResult);

    // Stampa le informazioni di TraceResult
    //printTraceResult(traceResult);

    exportTraceResult("TraceOutput.txt", traceResult);

    return 0;
}
