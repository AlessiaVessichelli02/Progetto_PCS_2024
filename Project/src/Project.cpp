#include "Project.hpp"
#include "DFN.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace DFNLibrary {

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
    istringstream convertN(line.substr(line.find(';')+1));
    convertN >> S;
    //leggo la terza riga
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
}

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
