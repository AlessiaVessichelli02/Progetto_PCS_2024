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

    Traces traces; // Inizializzazione della struttura Traces

    for (size_t i = 0; i < polygons.size(); ++i) {
        for (size_t j = i + 1; j < polygons.size(); ++j) {
            Plane plane1 = calculatePlaneEquation(polygons[i]);
            Plane plane2 = calculatePlaneEquation(polygons[j]);

            IntersectionLine intersectionLine = calculateIntersectionLine(plane1, plane2);
            cout << "Linea di intersezione tra i piani " << i << " e " << j << ": Punto (" << intersectionLine.point.x() << ", " << intersectionLine.point.y() << ", " << intersectionLine.point.z()
                 << ") Direzione (" << intersectionLine.direction.x() << ", " << intersectionLine.direction.y() << ", " << intersectionLine.direction.z() << ")" << endl << endl;

            calculateAndPrintIntersections(polygons, intersectionLine, i, j, traces); // Passaggio della struttura Traces per memorizzare le informazioni
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


    //checkTracePoints(traces, polygons);
    vector<TraceResult> results = checkTracePoints(traces, polygons);

    exportTraceResults(results, "TraceOutput.txt");

    return 0;
}
