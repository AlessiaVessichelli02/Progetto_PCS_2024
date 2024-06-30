//Questo file contiene le dichiarazioni dei test e le inclusioni necessarie.
//Includerà le librerie di Google Test, ifile di intestazione e le dichiarazioni delle funzioni di test.

#ifndef __TESTDFN_HPP
#define __TESTDFN_HPP

/*
#include "DFN.hpp"
#include "Project.hpp"
*/

#include <gtest/gtest.h>
#include "TestDFN_Utilities.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <math.h>
#include <list>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Eigen;

// Dichiarazioni delle funzioni di test

//Test per calculateDistance
TEST(CalculateDistanceTest, DistanceCalculation) {
    Point p1(0.0, 0.0, 0.0);
    Point p2(1.0, 1.0, 1.0);

    double expected_distance = sqrt(3.0);
    double actual_distance = calculateDistance(p1, p2);

    EXPECT_DOUBLE_EQ(actual_distance, expected_distance);
}

//Test per isPointOnEdge
TEST(IsPointOnEdgeTest, PointOnEdge) {
    Point p1(0.0, 0.0, 0.0);
    Point p2(2.0, 2.0, 0.0);
    Point p(1.0, 1.0, 0.0);

    bool result = isPointOnEdge(p1, p2, p);

    EXPECT_TRUE(result);
}

//Test per calculateLineIntersection
TEST(CalculateIntersectionBetweenLinesTest, ParallelLines) {
    Point p1(0.0, 0.0, 0.0);
    Vector3d direction1(1.0, 1.0, 0.0);

    Point p2(1.0, 1.0, 0.0);
    Vector3d direction2(1.0, 1.0, 0.0);

    try {
        Vector3d intersection = calculateIntersectionBetweenLines(p1, direction1, p2, direction2);

        // Se arriviamo qui, la funzione non ha rilevato che le linee sono parallele
        FAIL() << "Era prevista un'eccezione a causa di linee parallele, ma non ne è stata generata alcuna.";
    } catch (const std::runtime_error& e) {
        EXPECT_STREQ(e.what(), "Le rette sono parallele e non si intersecano.");
    } catch (const std::exception& e) {
        FAIL() << "Eccezione imprevista generata: " << e.what();
    }
}

//Test per doSegmentsOverlap
TEST(DoSegmentsOverlapTest, OverlapDetection) {
    Vector3d A(0.0, 0.0, 0.0);
    Vector3d B(2.0, 0.0, 0.0);
    Vector3d C(1.0, 1.0, 0.0);
    Vector3d D(1.0, -1.0, 0.0);

    Vector3d overlapStart, overlapEnd;
    bool result = doSegmentsOverlap(A, B, C, D, overlapStart, overlapEnd);

    EXPECT_TRUE(result);
}


#endif // TEST_DFN_HPP
