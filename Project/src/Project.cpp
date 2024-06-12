#include "Project.hpp"
#include "DFN.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

//namespace DFNLibrary {

void DFN::readDFN(const string& filename)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    string line;
    getline(file, line); // Read "# Number of Fractures"
    int numFractures;
    file >> numFractures;
    getline(file, line); // Read the remaining part of the line

    for (int i = 0; i < numFractures; ++i) {
        getline(file, line); // Read "# FracturesId; NumVertices"
        int id;
        int numVertices;
        char semicolon;
        file >> id >> semicolon >> numVertices;
        getline(file, line); // Read the remaining part of the line
        getline(file, line); // Read "# Vertices"

        Fracture fracture(id);

        // Temporary storage for vertices coordinates
        vector<double> xCoords(numVertices);
        vector<double> yCoords(numVertices);
        vector<double> zCoords(numVertices);

        // Read x coordinates
        for (int j = 0; j < numVertices; ++j) {
            file >> xCoords[j];
            if(j < numVertices - 1) file.ignore(1, ';'); // Ignore the semicolon
        }
        getline(file, line);

        // Read y coordinates
        for (int j = 0; j < numVertices; ++j) {
            file >> yCoords[j];
            if(j < numVertices - 1) file.ignore(1, ';'); // Ignore the semicolon
        }
        getline(file, line);

        // Read z coordinates
        for (int j = 0; j < numVertices; ++j) {
            file >> zCoords[j];
            if(j < numVertices - 1) file.ignore(1, ';'); // Ignore the semicolon
        }
        getline(file, line);

        // Combine coordinates into vertices
        for (int j = 0; j < numVertices; ++j) {
            Vertex vertex(xCoords[j], yCoords[j], zCoords[j]);
            fracture.addVertex(vertex);
        }

        fractures.push_back(fracture);

    }

    file.close();
}

void DFN::print() const {
    for(const auto& fracture : fractures) {
        fracture.print();
    }
}


//}

/*
bool ImportDFN(const string& inputFilePath)
{
    // Apertura del File
    ifstream file;
    file.open(inputFilePath);

    if (file.fail())
    {
        cerr<< "Errore nell'apertura del file"<< endl;
        return false;
    }

    //leggo le righe
    string line;

    //salto la prima riga
    while (!file.eof())
    {
        getline(file, line);

        // Skip Comment Line
        if(line[0] != '# Number of Fractures')
            break;
    }

    //leggo la seconda riga
    getline(file, line);
    istringstream convertN(line);
    convertN >> NumFractures;

    for(double i = 0; i < NumFractures; ++i){
        //salto la terza riga
        /*while (!file.eof())
        {
            getline(file, line);

            // Skip Comment Line
            if(line[0] != '# FractureId; NumVertices')
                break;
        }

        string header;
        getline(file,header);

        //leggo la quarta riga
        string line;
        getline(file, line);
        convertN.clear();
        convertN.str(line);
        convertN >> FractureID;
        convertN.ignore(1, ';');
        convertN >> NumVertices;

        //salto la quinta riga
        header.clear();
        getline(file,header);

        //leggo la sesto riga
        getline(file, line);
        convertN.clear();
        string line;
        getline(file, line);
        convertN.clear();
        convertN.str(line);
        convertN >> FractureID;
        convertN.ignore(1, ';');
        convertN >> NumVertices;
    }
*/
    /*
    istringstream convertN(line.substr(line.find(';')+1));
    getline(file, line);
    convertN.clear();
    convertN.str(line.substr(line.find(';')+1));
    convertN >> n;
    //salto la terza riga
    while (!file.eof())
    {
        getline(file, line);

        // Skip Comment Line
        if(line[0] != 'w;r')
            break;
    }
    */
//}

/*bool ImportMesh(const string &filepath, PolygonalMesh& mesh)
{
    //Verifico l'importazione di Cell0D e se è andata a buon fine stampo i marker
    if(!ImportCell0Ds(filepath + "/Cell0Ds.csv", mesh))
    {
        return false;
    }
    else
    {
        cout << "Cell0D marker:" << endl;
        if(mesh.Cell0DMarkers.size() == 0)
        {
            cout << "Non ci sono markers diversi da 0 per Cell0D" << endl;
        }
        else
        {
            for(auto it = mesh.Cell0DMarkers.begin(); it != mesh.Cell0DMarkers.end(); it++) //auto è un tipo
            {
                cout << "key:\t" << it -> first << "\t values:";
                for(const unsigned int id : it -> second)
                    cout << "\t" << id;
                cout << endl;
            }

        }

    }

    //Verifico l'importazione di Cell1D e se è andata a buon fine stampo i marker
    if(!ImportCell1Ds(filepath + "/Cell1Ds.csv", mesh))
    {
        return false;
    }
    else
    {
        cout << "Cell1D marker:" << endl;
        if(mesh.Cell1DMarkers.size() == 0)
        {
            cout << "Non ci sono markers diversi da 0 for Cell1D" << endl;
        }
        else
        {
            for(auto it = mesh.Cell1DMarkers.begin(); it != mesh.Cell1DMarkers.end(); it++)
            {
                cout << "key:\t" << it -> first << "\t values:";
                for(const unsigned int id : it -> second)
                    cout << "\t" << id;

                cout << endl;
            }
        }

    }*/

/*
#include "Project.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

namespace DFN {

double calculateLength(const Point3D& p1, const Point3D& p2) {
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2) + pow(p2.z - p1.z, 2));
}

vector<Fracture> readDFNFromFile(const string& filename) {
    ifstream infile(filename);
    string line;
    vector<Fracture> fractures;

    if (infile.is_open()) {
        // Leggi il numero di fratture
        getline(infile, line); // Ignora la linea di commento
        getline(infile, line);
        int numFractures = stoi(line);

        for (int i = 0; i < numFractures; ++i) {
            // Leggi l'identificativo della frattura e il numero dei vertici
            getline(infile, line); // Ignora la linea di commento
            getline(infile, line);
            istringstream iss(line);
            int fractureId, numVertices;
            char semicolon;
            iss >> fractureId >> semicolon >> numVertices;

            // Leggi i vertici
            getline(infile, line); // Ignora la linea di commento
            vector<Point3D> vertices(numVertices);
            for (int j = 0; j < numVertices; ++j) {
                infile >> vertices[j].x1 >> semicolon >> vertices[j].x2 >> semicolon >> vertices[j].x3 >> semicolon >> vertices[j].x4;
            }
            infile.ignore(numeric_limits<streamsize>::max(), '\n');

            fractures.push_back({fractureId, vertices});
        }
        infile.close();
    }

    return fractures;
}

void writeTracesToFile(const string& filename, const vector<Trace>& traces) {
    ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << "# Number of Traces\n" << traces.size() << "\n";
        outfile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2\n";
        for (const auto& trace : traces) {
            outfile << trace.id << "; " << trace.fractureId1 << "; " << trace.fractureId2 << "; "
                    << trace.point1.x << "; " << trace.point1.y << "; " << trace.point1.z << "; "
                    << trace.point2.x << "; " << trace.point2.y << "; " << trace.point2.z << "\n";
        }
        outfile.close();
    }
}

void writeFractureTracesToFile(const string& filename, const vector<vector<Trace>>& fractureTraces) {
    ofstream outfile(filename);
    if (outfile.is_open()) {
        for (size_t i = 0; i < fractureTraces.size(); ++i) {
            if (!fractureTraces[i].empty()) {
                outfile << "# FractureId; NumTraces\n" << fractureTraces[i][0].fractureId1 << "; " << fractureTraces[i].size() << "\n";
                outfile << "# TraceId; Tips; Length\n";
                for (const auto& trace : fractureTraces[i]) {
                    outfile << trace.id << "; " << trace.tips << "; " << trace.length << "\n";
                }
            }
        }
        outfile.close();
    }
}

bool doLinesIntersect(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4, Point3D& intersection) {
    // Calcola il determinante per determinare se le linee sono parallele
    double denom = (p4.z - p3.z) * (p2.y - p1.y) - (p4.y - p3.y) * (p2.z - p1.z);
    if (fabs(denom) < 1e-10) return false; // Le linee sono parallele

    double ua = ((p4.y - p3.y) * (p1.z - p3.z) - (p4.z - p3.z) * (p1.y - p3.y)) / denom;
    double ub = ((p2.y - p1.y) * (p1.z - p3.z) - (p2.z - p1.z) * (p1.y - p3.y)) / denom;

    // Calcola il punto di intersezione
    intersection.x = p1.x + ua * (p2.x - p1.x);
    intersection.y = p1.y + ua * (p2.y - p1.y);
    intersection.z = p1.z + ua * (p2.z - p1.z);

    return (ua >= 0.0 && ua <= 1.0 && ub >= 0.0 && ub <= 1.0);
}

bool isPointOnLineSegment(const Point3D& p, const Point3D& p1, const Point3D& p2) {
    double length = calculateLength(p1, p2);
    double d1 = calculateLength(p, p1);
    double d2 = calculateLength(p, p2);
    return fabs(length - (d1 + d2)) < 1e-10;
}

MyDFN::MyDFN(const string& filename) {
    fractures = readDFNFromFile(filename);
}

void MyDFN::calculateTraces() {
    // Implementazione del calcolo delle tracce
}

void MyDFN::findTracesBetweenFractures(const Fracture& f1, const Fracture& f2) {
    // Implementazione della ricerca delle tracce tra fratture
}

void MyDFN::classifyTrace(const Trace& trace, vector<Trace>& passTraces, vector<Trace>& nonPassTraces) {
    // Implementazione della classificazione delle tracce
}

void MyDFN::classifyAndSortTraces() {
    // Implementazione della classificazione e ordinamento delle tracce
}

void MyDFN::writeResults(const string& traceFilename, const string& fractureTraceFilename) {
    writeTracesToFile(traceFilename, traces);
    writeFractureTracesToFile(fractureTraceFilename, fractureTraces);
}

} // namespace DFN
*/
