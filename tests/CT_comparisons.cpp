#include "SimpleBounce/SimpleBounce.hpp"
#include "TestModels.hpp"
#include "catch.hpp"

TEST_CASE("Model 1") {
    Model1 model;
    simplebounce::BounceCalculator<1> bounce(&model, 3);

    std::array<double, 1> phiTV{5.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(52.366607861));

    BENCHMARK("solve Model1") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 1a") {
    Model1a model;
    simplebounce::BounceCalculator<1> bounce(&model, 3);

    std::array<double, 1> phiTV{1.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(1083.28));

    BENCHMARK("solve Model 1a") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 1b") {
    Model1b model;
    simplebounce::BounceCalculator<1> bounce(&model, 3);

    std::array<double, 1> phiTV{1.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(6.62391));

    BENCHMARK("solve Model 1b") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2") {
    Model2 model;
    simplebounce::BounceCalculator<2> bounce(&model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(20.82));

    BENCHMARK("solve Model 2") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2a") {
    Model2a model;
    simplebounce::BounceCalculator<2> bounce(&model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(1748.07));

    BENCHMARK("solve Model 2a") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2b") {
    Model2b model;
    simplebounce::BounceCalculator<2> bounce(&model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(4.44616));

    BENCHMARK("solve Model 2b") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 3") {
    Model3 model;
    simplebounce::BounceCalculator<3> bounce(&model, 3);
    std::array<double, 3> phiTV{1., 1., 1.}; // a point at which V<0
    std::array<double, 3> phiFV{0., 0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(21.9159));

    BENCHMARK("solve Model 3") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 4") {
    Model4 model;
    simplebounce::BounceCalculator<4> bounce(&model, 3);
    std::array<double, 4> phiTV{1., 1., 1., 1.}; // a point at which V<0
    std::array<double, 4> phiFV{0., 0., 0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(55.8148));

    BENCHMARK("solve Model 4") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 5") {
    Model5 model;
    simplebounce::BounceCalculator<5> bounce(&model, 3);
    std::array<double, 5> phiTV{1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 5> phiFV{0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(16.2475));

    BENCHMARK("solve Model 5") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 6") {
    Model6 model;
    simplebounce::BounceCalculator<6> bounce(&model, 3);
    std::array<double, 6> phiTV{1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 6> phiFV{0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(24.4173));

    BENCHMARK("solve Model 6") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 7") {
    Model7 model;
    simplebounce::BounceCalculator<7> bounce(&model, 3);

    std::array<double, 7> phiTV{1, 1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 7> phiFV{0, 0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(36.5977));

    BENCHMARK("solve Model 7") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 8") {
    Model8 model;
    simplebounce::BounceCalculator<8> bounce(&model, 3);

    std::array<double, 8> phiTV{1, 1, 1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 8> phiFV{0, 0, 0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == 0);
    REQUIRE(bounce.action() == Approx(45.9042));

    BENCHMARK("solve Model 8") { return bounce.solve(phiFV, phiTV); };
}