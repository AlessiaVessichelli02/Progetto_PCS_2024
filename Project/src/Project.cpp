#include "Project.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <map>
#include <iomanip>
#include <algorithm>
#include <tuple>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

bool readDFN(const string& filename, Fracture& frattura)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return 1;
    }

    file.precision(16);

    string line;

    // Leggo il numero di fratture
    getline(file, line); // Salto la prima riga 'NumFractures'
    getline(file, line); // Leggo il numero di fratture
    unsigned int numFractures;
    numFractures = stoi(line); // stoi() converte string in int

    frattura.NumFractures = numFractures; // Memorizzo il numero di fratture nell'apposita struttura

    //cout << "Numero di fratture: " << numFractures << endl;

    for (int i = 0; i < numFractures; ++i) {

        // Leggo l'ID delle fratture e il numero dei vertici
        getline(file, line); // Salto la linea '# FractureId; NumVertices'
        getline(file, line); // Leggo il numero di fratture
        stringstream ss(line);
        //char ch;
        int fractureId;
        int numVertices;

        ss >> fractureId;
        frattura.FractureID.push_back(fractureId);

        ss.ignore(numeric_limits<streamsize>::max(), ';'); // Ignora tutto ciò che c'è fino a ';' incluso
        ss >> numVertices; // Legge ciò che c'è dopo il ';'
        frattura.NumVertices.push_back(numVertices);

        //cout << "FractureId: " << fractureId << ", NumVertici: " << numVertices << endl;

        // Leggo i vertici
        getline(file, line); // Salto la linea '# Vertices'

        MatrixXd vertici(3, numVertices);

        for (int j = 0; j < 3; ++j) {
            getline(file, line);
            stringstream ssVertices(line);
            ssVertices.precision(16);
            vector<double> vertexRow; // Inizializzo un vettore per memorizzare i vertici riga per riga

            for (int k = 0; k < numVertices; ++k) {
                double value;
                ssVertices >> value;
                ssVertices.ignore(numeric_limits<streamsize>::max(), ';');
                vertexRow.push_back(value);
            }

            for (int l = 0; l < vertici.rows(); ++l) {
                bool isEmpty = true;  // Assume che la riga sia vuota

                // Verifica se la riga è effettivamente vuota
                for (int m = 0; m < vertici.cols(); ++m) {
                    if (vertici(l, m) != 0) {
                        isEmpty = false;
                        break;
                    }
                }

                // Se la riga è vuota, assegna i valori
                if (isEmpty) {
                    for (int m = 0; m < vertexRow.size(); ++m) {
                        vertici(l, m) = vertexRow[m];
                    }
                    break; // Esce dal ciclo una volta che la riga è stat riempita
                }

            }

        }

        frattura.Vertices.push_back(vertici);

        //cout << "Vertici: " << endl << vertici << endl;
    }

    file.close();
    return true;
}

/*
bool printFracture(const Fracture& frattura)
{
    cout << "Numero totale di fratture: " << frattura.NumFractures << endl;


    cout << "FractureID: ";
    for (unsigned int id : frattura.FractureID) {
        cout << id << " ";
    }
    cout << endl;

    cout << "NumVertices: ";
    for (unsigned int num : frattura.NumVertices) {
        cout << num << " ";
    }
    cout << endl;

    cout << "Vertices: " << endl;
    for (const MatrixXd& vertices : frattura.Vertices) {
        cout << vertices << endl;
    }

    for (unsigned int i = 0; i < frattura.NumFractures; ++i) {
        cout << "FractureID [" << i <<"]: " << frattura.FractureID[i] << endl;
        cout << "NumVertices [" << i <<"]: " << frattura.NumVertices[i] << endl;
        cout << "Vertices [" << i <<"]: " << endl << frattura.Vertices[i] << endl;
    }
    return true;
}
*/

//funzione che crea i quadrilateri utilizzando i vertici letti dal file
vector<Polygon> createPolygons(const Fracture& frattura) {
    vector<Polygon> polygons;
    for (const auto& vertices : frattura.Vertices) {
        Polygon poly;
        for (int i = 0; i < vertices.cols(); ++i) {
            Point point(vertices(0, i), vertices(1, i), vertices(2, i));
            poly.vertices.push_back(point);
        }
        polygons.push_back(poly);
    }
    return polygons;
}

// Funzione per calcolare l'equazione del piano
Plane calculatePlaneEquation(const Polygon& poly) {
    if (poly.vertices.size() < 3) {
        throw invalid_argument("Il Poligono deve avere almeno 3 vertici.");
    }

    // Prendi tre punti del poligono
    Point p1 = poly.vertices[0];
    Point p2 = poly.vertices[1];
    Point p3 = poly.vertices[2];

    // Calcola i vettori dei lati del poligono
    Point v1 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
    Point v2 = {p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};

    // Calcola il vettore normale al piano usando il prodotto vettoriale
    Point normal = {v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x};

    // Normalizza il vettore normale
    double length = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    normal.x /= length;
    normal.y /= length;
    normal.z /= length;

    // Calcola il coefficiente d usando uno dei punti del quadrilatero
    double d = -(normal.x * p1.x + normal.y * p1.y + normal.z * p1.z);

    return {normal.x, normal.y, normal.z, d};
}

//funzione che calcola la direzione della linea di intersezione tra due piani
Vector3d calculateIntersectionDirection(const Plane& plane1, const Plane& plane2)
{
    Vector3d normal1(plane1.a, plane1.b, plane1.c);
    Vector3d normal2(plane2.a, plane2.b, plane2.c);
    Vector3d direction = normal1.cross(normal2); //prodotto vettoriale
    return direction.normalized(); //normalizzazione
}

//funzione che trova un punto sulla linea di intersezione
Vector3d calculateIntersectionPoint(const Plane& plane1, const Plane& plane2, const Vector3d& direction)
{
    Matrix3d A;
    Vector3d b;

    // Popolo la matrice A e il vettore b
    A << plane1.a, plane1.b, plane1.c,
        plane2.a, plane2.b, plane2.c,
        direction.x(), direction.y(), direction.z();

    b << -plane1.d, -plane2.d, 0.0; //il prodotto scalare con la direzione è 0

    // Risolvo il sistema di equazioni
    Vector3d intersectionPoint = A.colPivHouseholderQr().solve(b); //effettuo la decomposizione QR con pivotazione su colonne della matrice A
    return intersectionPoint;
}

IntersectionLine calculateIntersectionLine(const Plane& plane1, const Plane& plane2)
{
    // Calcola la direzione della retta di intersezione
    Vector3d direction = calculateIntersectionDirection(plane1, plane2);

    // Calcola un punto sulla retta di intersezione
    Vector3d point = calculateIntersectionPoint(plane1, plane2, direction);

    return {direction, point};
}


// Funzione per calcolare l'equazione della retta passante per due punti
string calculateLineEquation(const Point& p1, const Point& p2)
{
    double a = p2.y - p1.y;
    double b = p1.x - p2.x;
    double c = a * p1.x + b * p1.y;
    return to_string(a) + "x + " + to_string(b) + "y = " + to_string(c);
}

VerticesLine extractLinesFromPolygons(const vector<Polygon>& polygons)
{
    VerticesLine verticesLine;
    for (const auto& poly : polygons) {
        vector<string> polyLines;
        for (size_t i = 0; i < poly.vertices.size(); ++i) {
            Point p1 = poly.vertices[i];
            Point p2 = poly.vertices[(i + 1) % poly.vertices.size()];
            polyLines.push_back(calculateLineEquation(p1, p2));
        }
        verticesLine.VerticesLines.push_back(polyLines);
    }
    return verticesLine;
}

Vector3d calculateLineIntersection(const Point& p1, const Vector3d& direction1, const Point& p2, const Vector3d& direction2) {
    Vector3d originVector(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    Vector3d crossProduct = direction1.cross(direction2);
    if (crossProduct.norm() < 1e-9) {
        throw runtime_error("Le rette sono parallele e non si intersecano.");
    }
    double t = originVector.cross(direction2).dot(crossProduct) / crossProduct.squaredNorm();
    double u = originVector.cross(direction1).dot(crossProduct) / crossProduct.squaredNorm();
    Vector3d intersectionPoint = Vector3d(p1.x, p1.y, p1.z) + direction1 * t;
    return intersectionPoint;
}

bool isPointOnSegment(const Point& p1, const Point& p2, const Vector3d& intersection) {
    double minX = min(p1.x, p2.x);
    double maxX = max(p1.x, p2.x);
    double minY = min(p1.y, p2.y);
    double maxY = max(p1.y, p2.y);
    double minZ = min(p1.z, p2.z);
    double maxZ = max(p1.z, p2.z);

    return (intersection.x() >= minX && intersection.x() <= maxX &&
            intersection.y() >= minY && intersection.y() <= maxY &&
            intersection.z() >= minZ && intersection.z() <= maxZ);
}

bool doSegmentsOverlap(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D, Vector3d& overlapStart, Vector3d& overlapEnd) {
    // Assicurati che A sia il punto iniziale e B il punto finale per il primo segmento
    Vector3d p1 = (A.x() <= B.x()) ? A : B;
    Vector3d q1 = (A.x() <= B.x()) ? B : A;

    // Assicurati che C sia il punto iniziale e D il punto finale per il secondo segmento
    Vector3d p2 = (C.x() <= D.x()) ? C : D;
    Vector3d q2 = (C.x() <= D.x()) ? D : C;

    // Verifica se c'è un'intersezione tra gli intervalli [p1, q1] e [p2, q2]
    if (std::max(p1.x(), p2.x()) <= std::min(q1.x(), q2.x())) {
        overlapStart = (p1.x() > p2.x()) ? p1 : p2;
        overlapEnd = (q1.x() < q2.x()) ? q1 : q2;
        return true;
    }

    return false;
}

void calculateAndPrintIntersections(const vector<Polygon>& polygons, const IntersectionLine& intersectionLine, size_t i, size_t j, Traces& traces)
{
    bool hasIntersectionI = false;
    bool hasIntersectionJ = false;

    vector<Vector3d> intersectionsI;
    vector<Vector3d> intersectionsJ;

    for (size_t l = 0; l < polygons[i].vertices.size(); ++l) {
        Point p1 = polygons[i].vertices[l];
        Point p2 = polygons[i].vertices[(l + 1) % polygons[i].vertices.size()];
        Vector3d direction = Vector3d(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z).normalized();

        try {
            Vector3d intersection = calculateLineIntersection(p1, direction, {intersectionLine.point.x(), intersectionLine.point.y(), intersectionLine.point.z()}, intersectionLine.direction);
            if (isPointOnSegment(p1, p2, intersection)) {
                hasIntersectionI = true;
                intersectionsI.push_back(intersection);
            }
        } catch (const exception& e) {
            cerr << "Errore nel calcolo dell'intersezione: " << e.what() << endl;
        }
    }

    for (size_t l = 0; l < polygons[j].vertices.size(); ++l) {
        Point p1 = polygons[j].vertices[l];
        Point p2 = polygons[j].vertices[(l + 1) % polygons[j].vertices.size()];
        Vector3d direction = Vector3d(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z).normalized();

        try {
            Vector3d intersection = calculateLineIntersection(p1, direction, {intersectionLine.point.x(), intersectionLine.point.y(), intersectionLine.point.z()}, intersectionLine.direction);
            if (isPointOnSegment(p1, p2, intersection)) {
                hasIntersectionJ = true;
                intersectionsJ.push_back(intersection);
            }
        } catch (const exception& e) {
            cerr << "Errore nel calcolo dell'intersezione: " << e.what() << endl;
        }
    }

    if (hasIntersectionI && hasIntersectionJ) {
        //vector<int> fractureIDs{i, j};
        vector<int> fractureIDs{static_cast<int>(i), static_cast<int>(j)};
        if (intersectionsI == intersectionsJ) {
            traces.traces[fractureIDs] = { {intersectionsI[0].x(), intersectionsI[0].y(), intersectionsI[0].z()}, {intersectionsI[1].x(), intersectionsI[1].y(), intersectionsI[1].z()} };
        } else {
            if (intersectionsI.size() >= 2 && intersectionsJ.size() >= 2) {
                Vector3d segmentI_start = intersectionsI[0];
                Vector3d segmentI_end = intersectionsI[1];
                Vector3d segmentJ_start = intersectionsJ[0];
                Vector3d segmentJ_end = intersectionsJ[1];

                Vector3d overlapStart, overlapEnd;
                if (doSegmentsOverlap(segmentI_start, segmentI_end, segmentJ_start, segmentJ_end, overlapStart, overlapEnd)) {
                    traces.traces[fractureIDs] = { {overlapStart.x(), overlapStart.y(), overlapStart.z()}, {overlapEnd.x(), overlapEnd.y(), overlapEnd.z()} };
                }
            }
        }
    }
}

void saveTracesToFile(const string& filename, const Traces& traces)
{
    ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    // Scrivi il numero di tracce nel file
    outputFile << "# Number of Traces" << endl;
    outputFile << traces.traces.size() << endl;

    int traceId; // Dichiarazione di traceId come intero

    // Ciclo per iterare attraverso le tracce nella mappa
    for (auto it = traces.traces.begin(); it != traces.traces.end(); ++it) {
        // Scrivi l'intestazione
        outputFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

        // Assegna l'ID della traccia
        traceId = distance(traces.traces.begin(), it);

        // Scrivi i dati della traccia
        outputFile << traceId << " " << it->first[0] << " " << it->first[1] << " ";

        // Scrivi i punti di intersezione
        for (const auto& point : it->second) {
            for (int i = 0; i < point.size(); ++i) {
                outputFile << point[i];
                if (i < point.size() - 1) {
                    outputFile << " ";
                }
            }
            outputFile << " ";
        }
        outputFile << endl;
    }

    outputFile.close();
}

// Funzione per verificare se un punto è su un segmento
bool isPointOnEdge(const Point& p1, const Point& p2, const Point& p)
{
    // Calcola il prodotto vettoriale per verificare se il punto è collineare con il segmento
    double crossProduct = (p.y - p1.y) * (p2.x - p1.x) - (p.x - p1.x) * (p2.y - p1.y);
    if (fabs(crossProduct) > 1e-6) return false;

    // Verifica se il punto è entro i limiti del segmento
    double minX = min(p1.x, p2.x);
    double maxX = max(p1.x, p2.x);
    double minY = min(p1.y, p2.y);
    double maxY = max(p1.y, p2.y);
    double minZ = min(p1.z, p2.z);
    double maxZ = max(p1.z, p2.z);

    return (p.x >= minX && p.x <= maxX &&
            p.y >= minY && p.y <= maxY &&
            p.z >= minZ && p.z <= maxZ);
}

// Funzione per calcolare la distanza euclidea tra due punti
double calculateDistance(const Point& p1, const Point& p2)
{
    return sqrt(pow(p2.x - p1.x, 2) +
                pow(p2.y - p1.y, 2) +
                pow(p2.z - p1.z, 2));
}

vector<TraceResult> checkTracePoints(const Traces& traces, const vector<Polygon>& polygons)
{
    vector<TraceResult> results;
    size_t traceId = 0;

    for (const auto& trace : traces.traces) {
        vector<int> fractureIDs = trace.first;
        vector<vector<double>> points = trace.second;

        Point p1(points[0][0], points[0][1], points[0][2]);
        Point p2(points[1][0], points[1][1], points[1][2]);
        double distance = calculateDistance(p1, p2);

        for (int fractureID : fractureIDs) {
            bool p1OnPolygon = false;
            bool p2OnPolygon = false;

            for (size_t i = 0; i < polygons[fractureID].vertices.size(); ++i) {
                Point polyPoint1 = polygons[fractureID].vertices[i];
                Point polyPoint2 = polygons[fractureID].vertices[(i + 1) % polygons[fractureID].vertices.size()];

                if (isPointOnEdge(polyPoint1, polyPoint2, p1)) {
                    p1OnPolygon = true;
                }
                if (isPointOnEdge(polyPoint1, polyPoint2, p2)) {
                    p2OnPolygon = true;
                }
            }

            bool isNonPassante = !(p1OnPolygon && p2OnPolygon);
            cout << "Per FractureId " << fractureID << " e TraceId " << traceId << " è "
                 << (isNonPassante ? "non passante" : "passante")
                 << " e la distanza fra i punti è " << distance << endl;

            results.push_back({static_cast<int>(fractureID), static_cast<int>(traceId), isNonPassante, distance});
        }

        traceId++;
    }

    return results;
}

// Funzione per confrontare due TraceResult in base alla lunghezza (distance)
bool compareByLength(const TraceResult& a, const TraceResult& b)
{
    return a.distance > b.distance; // Ordine decrescente
}

void exportTraceResults(const vector<TraceResult>& results, const string& filename) {
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Errore nell'apertura del file di output." << endl;
        return;
    }

    // Raggruppa i risultati per FractureId
    map<int, vector<TraceResult>> resultsByFracture;
    for (const auto& result : results) {
        resultsByFracture[result.fractureId].push_back(result);
    }

    // Scrivi i risultati nel file con il nuovo formato
    for (auto& [fractureId, traceResults] : resultsByFracture) {
        // Ordina i risultati per isNonPassante e poi per lunghezza (distance) in ordine decrescente
        sort(traceResults.begin(), traceResults.end(), [](const TraceResult& a, const TraceResult& b) {
            if (a.isNonPassante != b.isNonPassante) {
                return !a.isNonPassante && b.isNonPassante;  // Passanti (isNonPassante = false) prima di non passanti (isNonPassante = true)
            }
            return b.distance < a.distance;  // Ordinamento decrescente per lunghezza
        });

        outFile << "# FractureId; NumTraces\n";
        outFile << fractureId << "; " << traceResults.size() << "\n";
        outFile << "# TraceId; Tips; Length\n";
        for (const auto& result : traceResults) {
            outFile << result.traceId << "; " << (result.isNonPassante ? 1 : 0) << "; " << result.distance << "\n";
        }
    }

    outFile.close();
}
