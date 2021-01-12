#define SIMPLEBOUNCE_VERBOSE
#include "SimpleBounce/SimpleBounce.hpp"
#include <iostream>

class MyModel {
  public:
    // number of field degrees of freedom
    static constexpr std::size_t nPhi = 1;

    // potential for scalar field(s)
    double vpot(const double *phi) const {
        return phi[0] * phi[0] / 2. - phi[0] * phi[0] * phi[0] / 3.;
    }
    // first derivative(s) of potential
    void calcDvdphi(const double *phi, double *dvdphi) const {
        dvdphi[0] = phi[0] - phi[0] * phi[0];
    }
};

int main() {
    MyModel model{};
    std::array<double, 1> phiTV{10.}; // a point at which V<0
    std::array<double, 1> phiFV{0.};  // false vacuum

    // calculate the bounce solution
    auto bounceSolution{simplebounce::solve(model, phiFV, phiTV, 4)};

    // show the Euclidean action
    std::cout << "S_E = " << bounceSolution.action() << std::endl;

    return 0;
}
