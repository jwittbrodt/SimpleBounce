#include "SimpleBounce/SimpleBounce.hpp"
#include "TestModels.hpp"
#include "catch.hpp"

TEST_CASE("Model 1") {
    Model1 model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);

    std::array<double, 1> phiTV{5.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(52.4).epsilon(1e-2));

    BENCHMARK("solve Model1") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 1a") {
    Model1a model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);

    std::array<double, 1> phiTV{1.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(1.09e3).epsilon(5e-2));

    BENCHMARK("solve Model 1a") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 1b") {
    Model1b model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);

    std::array<double, 1> phiTV{1.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(6.63).epsilon(1e-2));

    BENCHMARK("solve Model 1b") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2") {
    Model2 model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(20.84).epsilon(1e-2));

    BENCHMARK("solve Model 2") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2a") {
    Model2a model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(1.75e3).epsilon(5e-2));

    BENCHMARK("solve Model 2a") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2b") {
    Model2b model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(4.456).epsilon(1e-2));

    BENCHMARK("solve Model 2b") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 3") {
    Model3 model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);
    std::array<double, 3> phiTV{1., 1., 1.}; // a point at which V<0
    std::array<double, 3> phiFV{0., 0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(21.95).epsilon(1e-2));

    BENCHMARK("solve Model 3") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 4") {
    Model4 model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);
    std::array<double, 4> phiTV{1., 1., 1., 1.}; // a point at which V<0
    std::array<double, 4> phiFV{0., 0., 0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(55.87).epsilon(1e-2));

    BENCHMARK("solve Model 4") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 5") {
    Model5 model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);
    std::array<double, 5> phiTV{1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 5> phiFV{0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(16.27).epsilon(1e-2));

    BENCHMARK("solve Model 5") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 6") {
    Model6 model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);
    std::array<double, 6> phiTV{1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 6> phiFV{0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(24.45).epsilon(1e-2));

    BENCHMARK("solve Model 6") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 7") {
    Model7 model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);

    std::array<double, 7> phiTV{1, 1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 7> phiFV{0, 0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(36.67).epsilon(1e-2));

    BENCHMARK("solve Model 7") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 8") {
    Model8 model;
    auto bounce = simplebounce::makeBounceCalculator(model, 3);

    std::array<double, 8> phiTV{1, 1, 1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 8> phiFV{0, 0, 0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(46.00).epsilon(1e-2));

    BENCHMARK("solve Model 8") { return bounce.solve(phiFV, phiTV); };
}
