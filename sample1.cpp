#include "SimpleBounce.hpp"
#include <iostream>
using namespace std;
using namespace simplebounce;

class MyModel : public GenericModel<1> {
  public:
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
    BounceCalculator<1> bounce(&model);
    bounce.verboseOn(); // verbose mode

    std::array<double, 1> phiTV{10.}; // a point at which V<0
    std::array<double, 1> phiFV{0.};  // false vacuum

    // calcualte the bounce solution
    bounce.solve(phiFV, phiTV);

    // show the results
    bounce.printBounce();

    // show the Euclidean action
    cout << "S_E = " << bounce.action() << endl;

    return 0;
}
