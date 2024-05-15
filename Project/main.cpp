#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "Project.hpp"
using namespace std;
int main()
{
    string inputFileName = "./FR3_data.txt";
    //leggo i vettori da file
    if (!ImportDFN()
    {
        cerr << "Errore durante l'importazione dei vettori" << endl;
        return -1;
    }
    else
        cout << "Importazione riuscita" << endl;
  return 0;
}
