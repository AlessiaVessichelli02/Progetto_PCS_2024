//Questo file contiene le dichiarazioni dei test e le inclusioni necessarie.
//Includerà le librerie di Google Test, ifile di intestazione e le dichiarazioni delle funzioni di test.

#ifndef __TESTDFN_HPP
#define __TESTDFN_HPP

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

//Test per isPointOnSegment
TEST(IsPointOnSegmentTest, PointOnSegment) {
    Point p1(0.0, 0.0, 0.0);
    Point p2(2.0, 2.0, 0.0);
    Vector3d v(1.0, 1.0, 0.0);

    bool result = isPointOnSegment(p1, p2, v);

    EXPECT_TRUE(result);
}

//Test per calculateIntersectionBetweenLines
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

// Test per calculateSphereRadius
TEST(CalculateSphereRadiusTest, RegularPolygonTest) {
    // Creiamo un poligono di esempio
    Polygon poly;
    poly.vertices = { {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}, {0.0, 2.0, 0.0} };

    // Calcoliamo il raggio della sfera
    double radius = calculateSphereRadius(poly);

    // Raggio atteso: in questo caso, la distanza dal baricentro (1.0, 1.0, 0.0) al vertice (2.0, 2.0, 0.0)
    double expectedRadius = sqrt(2.0); // Raggio atteso è la distanza tra (1.0, 1.0, 0.0) e (2.0, 2.0, 0.0)

    double tolDefault = 10 * numeric_limits<double>::epsilon();
    // Verifica con tolleranza
    EXPECT_NEAR(radius, expectedRadius, tolDefault);
}

// Test per doPolygonsIntersect
TEST(DoPolygonsIntersectTest, IntersectingPolygons) {
    Polygon poly1 = { {{0, 0, 0}, {2, 0, 0}, {1, 1.73, 0}} };//triangolo sul piano XY
    Polygon poly2 = { {{1, 0, 0}, {3, 0, 0}, {2, 1.73, 0}} };//triangolo traslato di (1,0,0)
    EXPECT_TRUE(doPolygonsIntersect(poly1, poly2));//i poligoni si intersecano
}

TEST(DoPolygonsIntersectTest, NonIntersectingPolygons) {
    Polygon poly1 = { {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}} };//triangolo sul piano XY
    Polygon poly2 = { {{2, 2, 0}, {3, 2, 0}, {2, 3, 0}} };//triangolo traslato di (2,2,0)
    EXPECT_FALSE(doPolygonsIntersect(poly1, poly2));//i poligoni non si intersecano
}

#endif // TEST_DFN_HPP

