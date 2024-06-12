#include "src/Project.hpp"
#include "src/DFN.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;
//using namespace DFNLibrary;

int main() {

    DFN dfn;
    dfn.readDFN("FR3_data.txt");
    dfn.print();

    return 0;
}


/*
int main()
{   
    string inputFileName = "FR3_data.txt";
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
*/

/*
#include "src/DFN.hpp"
#include "src/Project.hpp"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace DFN;

int main() {
    string inputFilename = "FR3_data.txt"; // Nome del file di input
    string traceOutputFilename = "traces.out";
    string fractureTraceOutputFilename = "fracture_traces.out";

    MyDFN dfn(inputFilename);
    dfn.calculateTraces();
    dfn.classifyAndSortTraces();
    dfn.writeResults(traceOutputFilename, fractureTraceOutputFilename);

    cout << "Tracce calcolate e risultati scritti nei file di output." << endl;

    return 0;
}
*/
