//Questo file contiene l'implementazione dei test e il punto di ingresso del programma di test.
//Definisce e implementa i casi di test e utilizza RUN_ALL_TESTS() per eseguire tutti i test definiti.

#include <gtest/gtest.h>
#include "src_test/TestDFN.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;
using namespace Eigen;

// Punto di ingresso per l'esecuzione dei test
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
