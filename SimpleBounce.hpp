#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace simplebounce {

template <class It>
double trapezoidalIntegrate(It beginValues, It endValues, double stepSize) {
    It curr = beginValues;
    It next = curr + 1;
    double result = 0.;
    while (next != endValues) {
        result += (*curr++ + *next++) * stepSize / 2.;
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////
template <std::size_t nPhi> class Scalarfield {
  private:
    int n_, dim_;
    double rmax_, dr_, drinv_;

    std::vector<std::array<double, nPhi>> phi_;
    std::vector<double> rinv_, r_dminusoneth_;

    void updateInfo() {
        for (int i = 1; i < n(); i++) {
            rinv_[i] = 1. / (dr_ * i);
            r_dminusoneth_[i] = std::pow(dr_ * i, dim_ - 1);
        }
    }

  public:
    Scalarfield(int n, double rmax, int dim)
        : n_{n}, dim_{dim}, rmax_{rmax}, dr_{rmax_ / (n_ - 1)}, drinv_{1 / dr_},
          phi_(n), rinv_(n_), r_dminusoneth_(n_) {
        updateInfo();
    }

    double phi(int i, int iphi) const { return phi_[i][iphi]; }

    double &phi(int i, int iphi) { return phi_[i][iphi]; }

    const std::array<double, nPhi> &phiVec(int i) const { return phi_[i]; }

    double lap(int i, int iphi) const {
        if (i == 0) {
            return 2. * (phi(1, iphi) - phi(0, iphi)) * drinv_ * drinv_ * dim_;
        } else {
            return (phi(i + 1, iphi) - 2. * phi(i, iphi) + phi(i - 1, iphi)) *
                       drinv_ * drinv_ +
                   (phi(i + 1, iphi) - phi(i - 1, iphi)) * 0.5 * drinv_ *
                       (dim_ - 1.) * rinv_[i];
        }
    }

    void setRmax(double rmax) {
        if (rmax > 0.) {
            rmax_ = rmax;
        } else {
            std::cerr << "!!! rmax should be positive value !!!" << std::endl;
            std::cerr << "!!! rmax is set to 1. !!!" << std::endl;
            rmax_ = 1.;
        }
        dr_ = rmax_ / (n_ - 1.);
        drinv_ = 1. / dr_;
        updateInfo();
    }

    void setDimension(int dim) {
        dim_ = dim;
        updateInfo();
    }

    void setN(int n) {
        n_ = n;
        dr_ = rmax_ / (n - 1.);
        drinv_ = 1. / dr_;
        phi_.resize(n_);
        rinv_.resize(n_);
        r_dminusoneth_.resize(n_);
        updateInfo();
    }

    int n() const { return n_; }
    int dim() const { return dim_; }
    double rmax() const { return rmax_; }
    double dr() const { return dr_; }
    double r_dminusoneth(const int i) const { return r_dminusoneth_[i]; }
};

template <std::size_t nphi, class VFunc, class DVFunc> class GenericModel {
    VFunc vFunc_;
    DVFunc dvFunc_;

  public:
    static constexpr std::size_t nPhi = nphi;

    GenericModel(VFunc vFunc, DVFunc dvFunc)
        : vFunc_{std::move(vFunc)}, dvFunc_{std::move(dvFunc)} {}

    double vPot(const std::array<double, nPhi> &phi) const {
        return vFunc_(phi);
    }
    std::array<double, nPhi>
    calcDvPot(const std::array<double, nPhi> &phi) const {
        return dvFunc_(phi);
    }
};

////////////////////////////////////////////////////////////////////////////////
template <class Model> class BounceCalculator {
  private:
    const Model &model_;

    Scalarfield<Model::nPhi> field_;

    double lambda;
    std::array<double, Model::nPhi> phiTV_;
    std::array<double, Model::nPhi> phiFV_;
    bool setVacuumDone;
    double VFV;
    bool verbose;

    // parameters for numerical calculation
    double safetyfactor;
    double maximumvariation;
    double xTV0;
    double width0;
    double derivMax;
    double tend0;
    double tend1;
    int maxN;

  public:
    static constexpr std::size_t nPhi = Model::nPhi;

    BounceCalculator(const Model &model) : field_{100, 1., 4}, model_{model} {

        // flags
        setVacuumDone = false;
        verbose = false;

        // parameters for numerical calculation
        safetyfactor = 0.9;
        maximumvariation = 0.01;
        xTV0 = 0.5;
        width0 = 0.05;
        derivMax = 1e-2;
        tend0 = 0.05;
        tend1 = 0.4;
        maxN = 1000;
    }

    double t() const {
        std::vector<double> integrand(field_.n() - 1);
        for (int i = 0; i < field_.n() - 1; i++) {
            integrand[i] = 0.;
            for (int iphi = 0; iphi < nPhi; iphi++) {
                integrand[i] += field_.r_dminusoneth(i) * -0.5 *
                                field_.phi(i, iphi) * field_.lap(i, iphi);
            }
        }
        return trapezoidalIntegrate(integrand.begin(), integrand.end(),
                                    field_.dr());
    }

    double v() const {
        std::vector<double> integrand(field_.n());
        for (int i = 0; i < field_.n(); i++) {
            integrand[i] =
                field_.r_dminusoneth(i) * (model_.vPot(field_.phiVec(i)) - VFV);
        }
        return trapezoidalIntegrate(integrand.begin(), integrand.end(),
                                    field_.dr());
    }

    double evolve(double dtau) {

        std::vector<std::array<double, nPhi>> laplacian(field_.n());
        // \nabla^2 \phi_iphi at r = r_i
        for (int i = 0; i < field_.n() - 1; i++) {
            for (int iphi = 0; iphi < nPhi; iphi++) {
                laplacian[i][iphi] = field_.lap(i, iphi);
            }
        }

        // integral1 : \int_0^\infty dr r^{d-1} \sum_i (\partial V /
        // \partial\phi_i) \nabla^2\phi_i integral2 : \int_0^\infty dr r^{d-1}
        // \sum_i (\partial V / \partial\phi_i)^2
        std::vector<double> integrand1(field_.n()), integrand2(field_.n());
        for (int i = 0; i < field_.n() - 1; i++) {
            integrand1[i] = 0.;
            integrand2[i] = 0.;
            auto dvdphi = model_.calcDvPot(field_.phiVec(i));
            for (int iphi = 0; iphi < nPhi; iphi++) {
                integrand1[i] +=
                    field_.r_dminusoneth(i) * dvdphi[iphi] * laplacian[i][iphi];
                integrand2[i] +=
                    field_.r_dminusoneth(i) * dvdphi[iphi] * dvdphi[iphi];
            }
        }
        integrand1[field_.n() - 1] = 0.;
        integrand2[field_.n() - 1] = 0.;

        // Eq. 9 of 1907.02417
        lambda = trapezoidalIntegrate(integrand1.begin(), integrand1.end(),
                                      field_.dr()) /
                 trapezoidalIntegrate(integrand2.begin(), integrand2.end(),
                                      field_.dr());

        // RHS of Eq. 8 of 1907.02417
        // phi at boundary is fixed to phiFV and will not be updated.
        std::vector<std::array<double, nPhi>> RHS(field_.n());
        for (int i = 0; i < field_.n() - 1; i++) {
            auto dvdphi = model_.calcDvPot(field_.phiVec(i));
            for (int iphi = 0; iphi < nPhi; iphi++) {
                RHS[i][iphi] = laplacian[i][iphi] - lambda * dvdphi[iphi];
            }
        }

        // if RHS of EOM at the origin is too big, smaller step is taken.
        double sum = 0.;
        for (int iphi = 0; iphi < nPhi; iphi++) {
            sum += RHS[0][iphi] * RHS[0][iphi];
        }
        double dtautilde = maximumvariation * fieldExcursion() / sqrt(sum);
        if (dtau < dtautilde) {
            dtautilde = dtau;
        }

        // flow by Eq. 8 of 1907.02417
        // phi at boundary is fixed to phiFV and will not be updated.
        for (int i = 0; i < field_.n() - 1; i++) {
            for (int iphi = 0; iphi < nPhi; iphi++) {
                field_.phi(i, iphi) += dtautilde * RHS[i][iphi];
            }
        }

        return lambda;
    }

    double tBounce() const {
        double area = field_.dim() * pow(M_PI, field_.dim() / 2.) /
                      tgamma(field_.dim() / 2. + 1.);
        return area * pow(lambda, field_.dim() / 2. - 1.) * t();
    }

    double vBounce() const {
        double area = field_.dim() * pow(M_PI, field_.dim() / 2.) /
                      tgamma(field_.dim() / 2. + 1.);
        return area * pow(lambda, field_.dim() / 2.) * v();
    }

    double action() const { return tBounce() + vBounce(); }

    double oneIfBounce() const {
        return (2. - field_.dim()) / field_.dim() * tBounce() / vBounce();
    }

    double rBounce(int i) const { return sqrt(lambda) * field_.dr() * i; }

    int setVacuum(const std::array<double, nPhi> &phiTV,
                  const std::array<double, nPhi> &phiFV) {
        if (model_.vPot(phiTV) > model_.vPot(phiFV)) {
            std::cerr
                << "!!! energy of true vacuum is larger than false vacuum !!!"
                << std::endl;
            return -1;
        }
        std::copy(phiTV.begin(), phiTV.end(), phiTV_.begin());
        std::copy(phiFV.begin(), phiFV.end(), phiFV_.begin());

        if (verbose) {
            std::cerr << "true and false vacua have been set." << std::endl;

            std::cerr << "\tfalse vacuum : ";
            std::cerr << "(";
            for (int iphi = 0; iphi < nPhi - 1; iphi++) {
                std::cerr << phiFV[iphi] << ", ";
            }
            std::cerr << phiFV[nPhi - 1] << ")\t";
            std::cerr << "V = " << model_.vPot(phiFV) << std::endl;

            std::cerr << "\ta point with smaller V : ";
            std::cerr << "(";
            for (int iphi = 0; iphi < nPhi - 1; iphi++) {
                std::cerr << phiTV[iphi] << ", ";
            }
            std::cerr << phiTV[nPhi - 1] << ")\t";
            std::cerr << "V = " << model_.vPot(phiTV) << std::endl;
        }

        VFV = model_.vPot(phiFV);
        setVacuumDone = true;

        return 0;
    }

    void setInitial(double frac, double width) {
        for (int i = 0; i < field_.n() - 1; i++) {
            for (int iphi = 0; iphi < nPhi; iphi++) {
                field_.phi(i, iphi) =
                    phiTV_[iphi] + (phiFV_[iphi] - phiTV_[iphi]) *
                                       (1. + tanh((i - field_.n() * frac) /
                                                  (field_.n() * width))) /
                                       2.;
            }
        }
        for (int iphi = 0; iphi < nPhi; iphi++) {
            field_.phi(field_.n() - 1, iphi) = phiFV_[iphi];
        }
        if (verbose) {
            std::cerr << "\t"
                      << "xTrueVacuum:\t" << frac << std::endl;
            std::cerr << "\t"
                      << "xWidth:\t" << width << std::endl;
            std::cerr << "\t"
                      << "V[phi] :\t" << v() << std::endl;
            std::cerr << "\t"
                      << "n :\t" << field_.n() << std::endl;
        }
    }

    double fieldExcursion() const {
        double normsquared = 0.;
        for (int iphi = 0; iphi < nPhi; iphi++) {
            normsquared +=
                pow(field_.phi(field_.n() - 1, iphi) - field_.phi(0, iphi), 2);
        }
        return sqrt(normsquared);
    }

    double derivativeAtBoundary() const {
        double normsquared = 0.;
        for (int iphi = 0; iphi < nPhi; iphi++) {
            normsquared += pow(field_.phi(field_.n() - 1, iphi) -
                                   field_.phi(field_.n() - 2, iphi),
                               2);
        }
        return sqrt(normsquared) / field_.dr();
    }

    int evolveUntil(double tauend) {

        // 1 + d + sqrt(1 + d) is maximum of absolute value of eigenvalue of
        // {{-2d, 2d},{ (1-d)/2 + 1, -2}}, which is discreitzed Laplacian for n
        // = 2. This value is 6 for d=3, and 7.23607 for d=4. The numerical
        // value of maximum of absolute value of eigenvalue of discretized
        // Laplacian for large n is 6 for d=3, and 7.21417 for d=4
        double dtau = 2. / (1. + field_.dim() + sqrt(1. + field_.dim())) *
                      pow(field_.dr(), 2) * safetyfactor;

        if (verbose) {
            std::cerr << "evolve until tau = " << tauend << ", (dtau = " << dtau
                      << ")" << std::endl;
        }

        for (double tau = 0.; tau < tauend; tau += dtau) {
            evolve(dtau);
            if (derivativeAtBoundary() * field_.rmax() / fieldExcursion() >
                derivMax) {
                return -1;
            }
        }
        return 0;
    }

    int solve() {
        if (!setVacuumDone) {
            std::cerr << "!!! correct vacua have not been set yet !!!"
                      << std::endl;
            return -1;
        }

        // make the bubble wall thin to get negative potential energy
        // if V is positive, make the wall thin.
        if (verbose) {
            std::cerr << "probing a thickness to get negative V[phi] ..."
                      << std::endl;
        }
        double xTV = xTV0;
        double width = width0;
        setInitial(xTV, width);
        while (v() > 0.) {
            width = width * 0.5;
            if (width * field_.n() < 1.) {
                if (verbose) {
                    std::cerr << "the current mesh is too sparse. increase the "
                                 "number of points."
                              << std::endl;
                }
                field_.setN(2 * field_.n());
            }
            if (field_.n() > maxN) {
                std::cerr << "!!! n became too large !!!" << std::endl;
                return -1;
            }
            setInitial(xTV, width);
        }

        // make the size of the bubble smaller enough than the size of the
        // sphere if dphi/dr at the boundary becomes too large during flow, take
        // smaller bounce configuration
        if (verbose) {
            std::cerr << "probing the size of the bounce configuration ..."
                      << std::endl;
        }
        while (evolveUntil(tend1 * std::pow(field_.rmax(), 2)) != 0) {
            if (verbose) {
                std::cerr << "the size of the bounce is too large. initial "
                             "condition is scale transformed."
                          << std::endl;
            }
            xTV = xTV * 0.5;
            width = width * 0.5;
            if (width * field_.n() < 1.) {
                if (verbose) {
                    std::cerr << "the current mesh is too sparse. increase the "
                                 "number of points."
                              << std::endl;
                }
                field_.setN(2 * field_.n());
            }
            // retry by using new initial condition
            setInitial(xTV, width);
            if (field_.n() > maxN) {
                std::cerr << "!!! n became too large !!!" << std::endl;
                return -1;
            }
        }

        return 0;
    }

    double getlambda() const { return lambda; }

    int printBounce() const {
        if (!setVacuumDone) {
            std::cerr << "!!! correct vacua have not been set yet !!!"
                      << std::endl;
            return -1;
        }

        std::cout << "# ";
        std::cout << "r\t";
        for (int iphi = 0; iphi < nPhi; iphi++) {
            std::cout << "phi[" << iphi << "]\t";
        }
        for (int iphi = 0; iphi < nPhi; iphi++) {
            std::cout << "RHS of EOM[" << iphi << "]\t";
        }
        std::cout << std::endl;

        for (int i = 0; i < field_.n(); i++) {
            std::cout << rBounce(i) << "\t";
            for (int iphi = 0; iphi < nPhi; iphi++) {
                std::cout << field_.phi(i, iphi) << "\t";
            }
            std::cout << std::endl;
        }

        return 0;
    }

    int printBounceDetails() const {
        if (!setVacuumDone) {
            std::cerr << "!!! correct vacua have not been set yet !!!"
                      << std::endl;
            return -1;
        }

        std::cout << "Bounce: " << std::endl;
        std::cout << "\tS =\t" << action() << std::endl;
        std::cout << "\tT =\t" << tBounce() << std::endl;
        std::cout << "\tV =\t" << vBounce() << std::endl;

        std::cout << "\tr = 0 : ";
        std::cout << "(";
        for (int iphi = 0; iphi < nPhi - 1; iphi++) {
            std::cout << field_.phi(0, iphi) << ", ";
        }
        std::cout << field_.phi(0, nPhi - 1) << ")\t";
        std::cout << "V = " << model_.vpot(field_.phiVec(0)) << std::endl;

        std::cout << "\tr = " << rBounce(field_.n() - 1) << " : ";
        std::cout << "(";
        for (int iphi = 0; iphi < nPhi - 1; iphi++) {
            std::cerr << field_.phi(field_.n() - 1, iphi) << ", ";
        }
        std::cout << field_.phi(field_.n() - 1, nPhi - 1) << ")\t";
        std::cout << "V = " << model_.vpot(field_.phiVec(field_.n() - 1))
                  << std::endl;

        std::cout << std::endl;

        std::cout << "Before rescaling: " << std::endl;
        std::cout << "\tT =\t" << t() << std::endl;
        std::cout << "\tV =\t" << v() << std::endl;
        std::cout << "\tlambda =\t" << lambda << std::endl;
        std::cout << std::endl;

        std::cout << "Goodness of the solution: " << std::endl;
        std::cout << "\tderiv\t"
                  << derivativeAtBoundary() / fieldExcursion() * field_.rmax()
                  << std::endl;
        std::cout << "\t(2-d)/d*T/V =\t" << oneIfBounce() << std::endl;

        return 0;
    }

    void setSafetyfactor(double x) { safetyfactor = x; }

    // set the parameter to determine the maximal variation of the field in
    // evolve()
    void setMaximumvariation(double x) { maximumvariation = x; }

    // set the initial value of the parameter to determine the initail
    // configuration
    void setXTV0(double x) { xTV0 = x; }

    // set the initial value of the parameter to determine the initail
    // configuration
    void setWidth0(double x) { width0 = x; }

    // set the maximal value of the derivative of field at the boundary
    void setDerivMax(double x) { derivMax = x; }

    // set tau0
    void setTend0(double x) { tend0 = x; }

    // set tau1
    void setTend1(double x) { tend1 = x; }

    // set the maximal value of the grid
    void setMaxN(int x) { maxN = x; }

    // turn on verbose mode
    void verboseOn() { verbose = true; }

    // turn off verbose mode
    void verboseOff() { verbose = false; }

    void setDimension(std::size_t dim) { field_.setDimension(dim); }
    void setN(std::size_t n) { field_.setN(n); }
};

} // namespace simplebounce
