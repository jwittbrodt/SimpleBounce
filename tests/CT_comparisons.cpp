#include "SimpleBounce/SimpleBounce.hpp"
#include "TestModels.hpp"
#include "catch.hpp"

TEST_CASE("Model 1") {
    Model1 model;
    simplebounce::BounceCalculator<1> bounce(&model, 3);

    std::array<double, 1> phiTV{5.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(52.3655920511));

    BENCHMARK("solve Model1") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 1a") {
    Model1a model;
    simplebounce::BounceCalculator<1> bounce(&model, 3);

    std::array<double, 1> phiTV{1.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(1083.1189124556));

    BENCHMARK("solve Model 1a") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 1b") {
    Model1b model;
    simplebounce::BounceCalculator<1> bounce(&model, 3);

    std::array<double, 1> phiTV{1.}; // a point at which V<0
    std::array<double, 1> phiFV{0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(6.6237816762));

    BENCHMARK("solve Model 1b") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2") {
    Model2 model;
    simplebounce::BounceCalculator<2> bounce(&model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(20.8197046624));

    BENCHMARK("solve Model 2") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2a") {
    Model2a model;
    simplebounce::BounceCalculator<2> bounce(&model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(1748.1421154732));

    BENCHMARK("solve Model 2a") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 2b") {
    Model2b model;
    simplebounce::BounceCalculator<2> bounce(&model, 3);
    std::array<double, 2> phiTV{1., 1.}; // a point at which V<0
    std::array<double, 2> phiFV{0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(4.4459478241));

    BENCHMARK("solve Model 2b") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 3") {
    Model3 model;
    simplebounce::BounceCalculator<3> bounce(&model, 3);
    std::array<double, 3> phiTV{1., 1., 1.}; // a point at which V<0
    std::array<double, 3> phiFV{0., 0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(21.9150686391));

    BENCHMARK("solve Model 3") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 4") {
    Model4 model;
    simplebounce::BounceCalculator<4> bounce(&model, 3);
    std::array<double, 4> phiTV{1., 1., 1., 1.}; // a point at which V<0
    std::array<double, 4> phiFV{0., 0., 0., 0.}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(55.8136776842));

    BENCHMARK("solve Model 4") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 5") {
    Model5 model;
    simplebounce::BounceCalculator<5> bounce(&model, 3);
    std::array<double, 5> phiTV{1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 5> phiFV{0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(16.2470608546));

    BENCHMARK("solve Model 5") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 6") {
    Model6 model;
    simplebounce::BounceCalculator<6> bounce(&model, 3);
    std::array<double, 6> phiTV{1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 6> phiFV{0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(24.4165889149));

    BENCHMARK("solve Model 6") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 7") {
    Model7 model;
    simplebounce::BounceCalculator<7> bounce(&model, 3);

    std::array<double, 7> phiTV{1, 1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 7> phiFV{0, 0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(36.5962442623));

    BENCHMARK("solve Model 7") { return bounce.solve(phiFV, phiTV); };
}

TEST_CASE("Model 8") {
    Model8 model;
    simplebounce::BounceCalculator<8> bounce(&model, 3);

    std::array<double, 8> phiTV{1, 1, 1, 1, 1, 1, 1, 1}; // a point at which V<0
    std::array<double, 8> phiFV{0, 0, 0, 0, 0, 0, 0, 0}; // false vacuum
    REQUIRE(bounce.solve(phiFV, phiTV) == Approx(45.9021687082));

    BENCHMARK("solve Model 8") { return bounce.solve(phiFV, phiTV); };
}
