#pragma once

#include <algorithm>
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
template <std::size_t nPhi> class GenericModel {
  public:
    virtual double vpot(const double *phi) const = 0;
    virtual void calcDvdphi(const double *phi, double *dvdphi) const = 0;
};

struct Parameters {
    // initial grid size
    std::size_t initialN = 100;
    // maximal grid size
    std::size_t maxN = 1000;
    // radius at the boundary
    double rMax = 1.;
    double safetyfactor = 0.9;
    // maximal variation of the field in evolve()
    double maximumvariation = 0.01;
    // initial value of the parameter to determine the initail configuration
    double xTV0 = 0.5;
    // initial value of the parameter to determine the initail configuration
    double width0 = 0.05;
    // maximal value of the derivative of the field at the boundary
    double derivMax = 1e-2;
    // set tau1
    double tend1 = 0.4;
};

template <std::size_t nPhi> class FieldConfiguration {
    std::size_t n_, dim_;
    double rMax_, dr_, drinv_;
    double *phi_;
    std::vector<double> r_dminusoneth_;

    // radius : r_i = i * dr
    double r(int i) const { return dr_ * i; }

    // calculate and save pow(r(i), dim-1)
    void updateInfo() {
        r_dminusoneth_.resize(n_);
        for (int i = 0; i < n(); i++) {
            r_dminusoneth_[i] = pow(r(i), dim_ - 1);
        }
    }

  public:
    FieldConfiguration(std::size_t n, std::size_t dim, double rMax)
        : n_{n}, dim_{dim}, rMax_{rMax}, dr_{rMax_ / (n_ - 1)},
          drinv_{1. / dr_}, phi_{new double[n_ * nPhi]}, r_dminusoneth_{} {
        r_dminusoneth_.reserve(n_);
        for (int i = 0; i < n_; i++) {
            r_dminusoneth_.emplace_back(std::pow(r(i), dim_ - 1));
        }
    }

    // set the initial configuration
    // See Eq. 11 in the manual.
    void setInitial(const double frac, const double width,
                    const std::array<double, nPhi> &phiFV,
                    const std::array<double, nPhi> &phiTV) {
        for (int i = 0; i < n_ - 1; i++) {
            for (int iphi = 0; iphi < nPhi; iphi++) {
                setPhi(i, iphi,
                       phiTV[iphi] +
                           (phiFV[iphi] - phiTV[iphi]) *
                               (1. + tanh((i - n_ * frac) / (n_ * width))) /
                               2.);
            }
        }
        for (int iphi = 0; iphi < nPhi; iphi++) {
            setPhi(n_ - 1, iphi, phiFV[iphi]);
        }
    }

    ~FieldConfiguration() { delete[] phi_; }

    // return the value of scalar field phi_iphi at r_i
    double phi(const int i, const int iphi) const {
        return phi_[i * nPhi + iphi];
    }

    // set the value of scalar field phi_iphi to phi
    void setPhi(const int i, const int iphi, const double phi) {
        phi_[i * nPhi + iphi] = phi;
    }

    // add phi to the value of scalar field phi_iphi
    void addToPhi(const int i, const int iphi, const double phi) {
        phi_[i * nPhi + iphi] += phi;
    }

    // return the address of phi_0 at r_i.
    // phi_0, phi_1, phi_2, ... can be obtained by x[0], x[1], x[2], ... after x
    // = phivec(i)
    double *phivec(const int i) const { return &phi_[i * nPhi + 0]; }

    // set the dimension of the Euclidean space
    void setDimension(const int dim) {
        dim_ = dim;
        updateInfo();
    }

    // set the number of grid. grid spacing dr is consistently changed.
    void setN(const int n) {
        n_ = n;
        dr_ = rMax_ / (n - 1.);
        drinv_ = 1. / dr_;
        delete[] phi_;
        phi_ = new double[n * nPhi];
        updateInfo();
    }

    // return the number of the grid
    int n() const { return n_; }

    // return the dimension of space
    int dim() const { return dim_; }

    // return the lattice spacing
    double dr() const { return dr_; }

    // return pow(r(i), dim-1)
    double r_dminusoneth(const int i) const { return r_dminusoneth_[i]; }

    // Laplacian in radial coordinate : \nabla^2 \phi = d^2 phi / dr^2 + (d-1)/r
    // * dphi/dr note that rinv_[i] = 1/r(i). See Eq. 9 in the manual
    double lap(const int i, const int iphi) const {
        if (i == 0) {
            return 2. * (phi(1, iphi) - phi(0, iphi)) * drinv_ * drinv_ * dim_;
        } else {
            return (phi(i + 1, iphi) - 2. * phi(i, iphi) + phi(i - 1, iphi)) *
                       drinv_ * drinv_ +
                   (phi(i + 1, iphi) - phi(i - 1, iphi)) * 0.5 * drinv_ *
                       (dim_ - 1.) / r(i);
        }
    }

    // field excursion from the origin to the infinity
    double fieldExcursion() const {
        double normsquared = 0.;
        for (int iphi = 0; iphi < nPhi; iphi++) {
            normsquared += pow(phi(n_ - 1, iphi) - phi(0, iphi), 2);
        }
        return sqrt(normsquared);
    }

    // derivative of scalar field at boundary
    double derivativeAtBoundary() const {
        double normsquared = 0.;
        for (int iphi = 0; iphi < nPhi; iphi++) {
            normsquared += pow(phi(n_ - 1, iphi) - phi(n_ - 2, iphi), 2);
        }
        return sqrt(normsquared) * drinv_;
    }
};

// kinetic energy of the configuration
// \int_0^\infty dr r^{d-1} \sum_i (-1/2) \phi_i \nabla^2\phi_i
template <std::size_t nPhi>
double kineticEnergy(const FieldConfiguration<nPhi> &field) {
    std::vector<double> integrand(field.n() - 1);
    for (int i = 0; i < field.n() - 1; i++) {
        for (int iphi = 0; iphi < nPhi; iphi++) {
            integrand[i] -= field.r_dminusoneth(i) * field.phi(i, iphi) *
                            field.lap(i, iphi) * 0.5;
        }
    }
    return trapezoidalIntegrate(integrand.begin(), integrand.end(), field.dr());
}

// potential energy of the configuration
// \int_0^\infty dr r^{d-1} V(\phi)
template <std::size_t nPhi, class Model>
double potentialEnergy(const FieldConfiguration<nPhi> &field,
                       const Model *model, double zeroPot) {
    std::vector<double> integrand(field.n());
    for (int i = 0; i < field.n(); i++) {
        integrand[i] =
            field.r_dminusoneth(i) * (model->vpot(field.phivec(i)) - zeroPot);
    }
    return trapezoidalIntegrate(integrand.begin(), integrand.end(), field.dr());
}

template <std::size_t nPhi> class BounceCalculator {
  private:
    Parameters params_;
    std::size_t dim_;

    double lambda;

    GenericModel<nPhi> *model_;
    double VFV;
    bool verbose = false;

  public:
    BounceCalculator(GenericModel<nPhi> *model, std::size_t dim = 4,
                     Parameters params = {})
        : params_{params}, dim_{dim}, model_{model} {}

    // main routine to get the bounce solution
    // See Fig. 1 of the manual
    double solve(const std::array<double, nPhi> &phiFV,
                 const std::array<double, nPhi> &phiTV) {
        if (model_->vpot(phiTV.data()) > model_->vpot(phiFV.data())) {
            std::cerr
                << "!!! energy of true vacuum is larger than false vacuum !!!"
                << std::endl;
            return -1;
        }
        VFV = model_->vpot(phiFV.data());

        // make the bubble wall thin to get negative potential energy
        // if V is positive, make the wall thin.
        if (verbose) {
            std::cerr << "probing a thickness to get negative V[phi] ..."
                      << std::endl;
        }
        double xTV = params_.xTV0;
        double width = params_.width0;
        FieldConfiguration<nPhi> field{params_.initialN, dim_, params_.rMax};
        field.setInitial(xTV, width, phiFV, phiTV);
        while (potentialEnergy(field, model_, VFV) > 0.) {
            width = width * 0.5;
            if (width * field.n() < 1.) {
                if (verbose) {
                    std::cerr << "the current mesh is too sparse. increase the "
                                 "number of points."
                              << std::endl;
                }
                field.setN(2 * field.n());
            }
            if (field.n() > params_.maxN) {
                std::cerr << "!!! n became too large !!!" << std::endl;
                return -1;
            }
            field.setInitial(xTV, width, phiFV, phiTV);
        }
        // make the size of the bubble smaller enough than the size of the
        // sphere if dphi/dr at the boundary becomes too large during flow,
        // take smaller bounce configuration
        if (verbose) {
            std::cerr << "probing the size of the bounce configuration ..."
                      << std::endl;
        }
        while (evolveUntil(field,
                           params_.tend1 * params_.rMax * params_.rMax) != 0) {
            if (verbose) {
                std::cerr << "the size of the bounce is too large. initial "
                             "condition is scale transformed."
                          << std::endl;
            }
            xTV = xTV * 0.5;
            width = width * 0.5;
            if (width * field.n() < 1.) {
                if (verbose) {
                    std::cerr << "the current mesh is too sparse. increase the "
                                 "number of points."
                              << std::endl;
                }
                field.setN(2 * field.n());
            }
            // retry by using new initial condition
            field.setInitial(xTV, width, phiFV, phiTV);
            if (field.n() > params_.maxN) {
                std::cerr << "!!! n became too large !!!" << std::endl;
                return -1;
            }
        }
        return action(field);
    }

    // evolve the configuration from tau = 0 to tau = tauend
    int evolveUntil(FieldConfiguration<nPhi> &field, const double tauend) {

        // 1 + d + sqrt(1 + d) is maximum of absolute value of eigenvalue of
        // {{-2d, 2d},{ (1-d)/2 + 1, -2}}, which is discreitzed Laplacian for n
        // = 2. This value is 6 for d=3, and 7.23607 for d=4. The numerical
        // value of maximum of absolute value of eigenvalue of discretized
        // Laplacian for large n is 6 for d=3, and 7.21417 for d=4
        double dtau = 2. / (1. + dim_ + sqrt(1. + dim_)) * pow(field.dr(), 2) *
                      params_.safetyfactor;

        if (verbose) {
            std::cerr << "evolve until tau = " << tauend << ", (dtau = " << dtau
                      << ")" << std::endl;
        }

        for (double tau = 0.; tau < tauend; tau += dtau) {
            evolve(field, dtau);
            if (field.derivativeAtBoundary() * params_.rMax /
                    field.fieldExcursion() >
                params_.derivMax) {
                return -1;
            }
        }
        return 0;
    }

    // evolve the configuration by dtau
    double evolve(FieldConfiguration<nPhi> &field, double dtau) {
        const auto n = field.n();
        double laplacian[n][nPhi];

        // \nabla^2 \phi_iphi at r = r_i
        for (int i = 0; i < n - 1; i++) {
            for (int iphi = 0; iphi < nPhi; iphi++) {
                laplacian[i][iphi] = field.lap(i, iphi);
            }
        }

        // integral1 : \int_0^\infty dr r^{d-1} \sum_i (\partial V /
        // \partial\phi_i) \nabla^2\phi_i integral2 : \int_0^\infty dr r^{d-1}
        // \sum_i (\partial V / \partial\phi_i)^2
        double integrand1[n], integrand2[n];
        for (int i = 0; i < n - 1; i++) {
            integrand1[i] = 0.;
            integrand2[i] = 0.;
            double dvdphi[nPhi];
            model_->calcDvdphi(field.phivec(i), dvdphi);
            for (int iphi = 0; iphi < nPhi; iphi++) {
                integrand1[i] +=
                    field.r_dminusoneth(i) * dvdphi[iphi] * laplacian[i][iphi];
                integrand2[i] +=
                    field.r_dminusoneth(i) * dvdphi[iphi] * dvdphi[iphi];
            }
        }
        integrand1[n - 1] = 0.;
        integrand2[n - 1] = 0.;

        // Eq. 9 of 1907.02417
        lambda = trapezoidalIntegrate(&integrand1[0], &integrand1[0] + n,
                                      field.dr()) /
                 trapezoidalIntegrate(&integrand2[0], &integrand2[0] + n,
                                      field.dr());

        // RHS of Eq. 8 of 1907.02417
        // phi at boundary is fixed to phiFV and will not be updated.
        double RHS[n][nPhi];
        for (int i = 0; i < n - 1; i++) {
            double dvdphi[nPhi];
            model_->calcDvdphi(field.phivec(i), dvdphi);
            for (int iphi = 0; iphi < nPhi; iphi++) {
                RHS[i][iphi] = laplacian[i][iphi] - lambda * dvdphi[iphi];
            }
        }

        // if RHS of EOM at the origin is too big, smaller step is taken.
        double sum = 0.;
        for (int iphi = 0; iphi < nPhi; iphi++) {
            sum += RHS[0][iphi] * RHS[0][iphi];
        }
        double dtautilde =
            params_.maximumvariation * field.fieldExcursion() / sqrt(sum);
        if (dtau < dtautilde) {
            dtautilde = dtau;
        }

        // flow by Eq. 8 of 1907.02417
        // phi at boundary is fixed to phiFV and will not be updated.
        for (int i = 0; i < n - 1; i++) {
            for (int iphi = 0; iphi < nPhi; iphi++) {
                field.addToPhi(i, iphi, dtautilde * RHS[i][iphi]);
            }
        }

        return lambda;
    }

    // Kinetic energy of the bounce
    // The bounce configuration can be obtained from Eq. 15 of 1907.02417
    double tBounce(const FieldConfiguration<nPhi> &field) const {
        double area = dim_ * pow(M_PI, dim_ / 2.) / tgamma(dim_ / 2. + 1.);
        return area * pow(lambda, dim_ / 2. - 1.) * kineticEnergy(field);
    }

    // Potential energy of the bounce
    // The bounce configuration can be obtained from Eq. 15 of 1907.02417
    double vBounce(const FieldConfiguration<nPhi> &field) const {
        double area = dim_ * pow(M_PI, dim_ / 2.) / tgamma(dim_ / 2. + 1.);
        return area * pow(lambda, dim_ / 2.) *
               potentialEnergy(field, model_, VFV);
    }

    // Euclidean action in d-dimensional space
    // The bounce configuration can be obtained from Eq. 15 of 1907.02417
    double action(const FieldConfiguration<nPhi> &field) const {
        return tBounce(field) + vBounce(field);
    }
    // turn on verbose mode
    void verboseOn() { verbose = true; }

    // turn off verbose mode
    void verboseOff() { verbose = false; }
};

} // namespace simplebounce
