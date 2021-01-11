#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>
#ifdef SIMPLEBOUNCE_VERBOSE
#include <iostream>
#endif

namespace simplebounce {

namespace Detail {

const auto squaredSum = [](double sum, double x) {
    return sum + std::pow(x, 2);
};

template <class It>
double trapezoidalIntegrate(It beginValues, It endValues, double stepSize) {
    double borderVal = (*beginValues++ + *endValues--) / 2.;
    return stepSize * std::accumulate(beginValues, endValues, borderVal);
}

//! Computes the euclidean vector norm \f$ | \vec{x1} - \vec{x2} | \f$.
template <std::size_t n>
double normDifference(const std::array<double, n> &x1,
                      const std::array<double, n> &x2) {
    return std::sqrt(std::inner_product(x1.begin(), x1.end(), x2.begin(), 0.,
                                        Detail::squaredSum,
                                        std::minus<double>{}));
}
} // namespace Detail

struct Parameters {
    // initial grid size
    std::size_t initialN = 100;
    // radius at the boundary
    double rMax = 1.;
    double safetyfactor = 0.9;
    // maximal variation of the field in evolve()
    double maximumvariation = 0.01;
    // initial value of the parameter InitBounceParams::r0Frac
    double r0Frac = 0.5;
    // initial value of the parameter InitBounceParams::width
    double width = 0.05;
    // maximal value of the derivative of the field at the boundary
    double derivMax = 1e-2;
    // set tau1
    double tend1 = 0.4;
};

template <std::size_t nPhi> class FieldConfiguration;

template <std::size_t nPhi> class InitialBounceConfiguration {
    std::array<double, nPhi> phiFV_;
    std::array<double, nPhi> phiTV_;
    double width_;
    double r0byR_;

    static constexpr double half = 1 / 2.;
    static constexpr double almost0 = 0.1 / FieldConfiguration<nPhi>::maxN;
    static constexpr double almost1 = 1 - almost0;

  public:
    InitialBounceConfiguration(const std::array<double, nPhi> &phiFV,
                               const std::array<double, nPhi> &phiTV,
                               double width, double r0byR)
        : phiFV_{phiFV}, phiTV_{phiTV}, width_{width}, r0byR_{r0byR} {}

    std::array<double, nPhi> phi(double rbyR) const noexcept {
        if (rbyR > almost1) {
            return phiFV_;
        }
        if (rbyR < almost0) {
            return phiTV_;
        }
        const auto interpolatingField =
            [scaleFac = (1. + std::tanh((rbyR - r0byR_) / width_)) *
                        half](double phiT, double phiF) {
                return phiT + scaleFac * (phiF - phiT);
            };
        std::array<double, nPhi> res;
        std::transform(phiTV_.begin(), phiTV_.end(), phiFV_.begin(),
                       res.begin(), interpolatingField);
        return res;
    }

    std::vector<std::array<double, nPhi>> onGrid(std::size_t n) const {
        auto phiGrid = std::vector<std::array<double, nPhi>>(n);
        onGrid(phiGrid);
        return phiGrid;
    }

    void onGrid(std::vector<std::array<double, nPhi>> &grid) const {
        const auto n = grid.size();
        grid.front() = phiTV_;
        const auto delta = 1. / static_cast<double>(n - 1);
        for (std::size_t i = 1; i != n - 1; ++i) {
            grid[i] = phi(i * delta);
        }
        grid.back() = phiFV_;
    }

    double width() const noexcept { return width_; }
    double &width() noexcept { return width_; }
    double &r0byR() noexcept { return r0byR_; }

    //! Determine if this initial configuration has a negative potential energy
    //! for this model and dimensionality.
    //!
    //! Uses a dedicated trapezoidal integration adapted from
    //! boost::math::quadrature::trapezoidal to perform the \f$r\f$ integral
    //! from 0 to 1. Stops the integration as soon as the absolute value is
    //! larger than margin * the integration error or after maxRefinements.
    //! Returns true if the resulting value is significantly less than
    //! maxVal.
    template <class Model>
    bool negativePotentialEnergy(const Model &model,
                                 std::size_t dim) const noexcept {
        static constexpr std::size_t maxRefinements = 10;
        static constexpr double margin = 5;
        static constexpr double maxVal = -1e-5;

        const double zeroPot = model.vpot(phiFV_.data());
        const auto integrand = [dim, zeroPot, model, this](double r) {
            return std::pow(r, dim - 1) * (model.vpot(phi(r).data()) - zeroPot);
        };

        // the integrand vanishes at the endpoints by construction
        auto I0 = 0.;
        auto h = half;
        auto I1 = integrand(h) * h;

        // The recursion is:
        // I_k = 1/2 I_{k-1} + 1/2^k \sum_{j=1; j odd, j < 2^k} f(a+j(b-a)/2^k)
        auto k = 2ul;
        auto error = std::abs(I0 - I1);
        while (k < 4 || (k < maxRefinements && margin * error > std::abs(I1))) {
            I0 = I1;
            I1 = I0 * half;
            const auto p = std::size_t{1} << k;
            h *= half;
            auto sum = 0.;
            for (std::size_t j = 1; j < p; j += 2) {
                auto y = integrand(j * h);
                sum += y;
            }
            I1 += sum * h;
            ++k;
            error = std::abs(I0 - I1);
        }
        return I1 < maxVal - error;
    }
};

template <std::size_t nPhi> class FieldConfiguration {
  public:
    // maximal grid size
    static constexpr std::size_t maxN = 1000;

  private:
    std::size_t n_, dim_;
    double rMax_;
    std::vector<std::array<double, nPhi>> phi_;
    std::array<double, nPhi> phiFV_;
    std::array<double, nPhi> phiTV_;

    double dr_ = rMax_ / (n_ - 1);
    double drinv_ = 1 / dr_;
    std::vector<double> rToDimMin1_;

    // calculate \f$ r_i^{dim-1} \f$
    std::vector<double> calcRToDimMin1() const {
        std::vector<double> result(n_, std::pow(dr_, dim_ - 1));
        for (std::size_t i = 0; i != n_; ++i) {
            result[i] *= std::pow(i, dim_ - 1);
        }
        return result;
    }

    std::size_t
    determineGridSize(const InitialBounceConfiguration<nPhi> &initBounce) {
        const auto n = std::max(
            std::size_t{100}, static_cast<std::size_t>(1 / initBounce.width()));
        if (n > maxN) {
            throw(std::runtime_error("Maximum grid size exceeded"));
        }
        return n;
    }

  public:
    std::array<double, nPhi> &operator[](std::size_t i) { return phi_[i]; }
    const std::array<double, nPhi> &operator[](std::size_t i) const {
        return phi_[i];
    }

    FieldConfiguration(const InitialBounceConfiguration<nPhi> &initBounce,
                       std::size_t dim, double rMax)
        : n_{determineGridSize(initBounce)}, dim_{dim}, rMax_{rMax},
          phi_{initBounce.onGrid(n_)}, phiFV_{initBounce.phi(1.)},
          phiTV_{initBounce.phi(0.)}, rToDimMin1_{calcRToDimMin1()} {}

    void
    resetInitialBounce(const InitialBounceConfiguration<nPhi> &initBounce) {
        auto n = determineGridSize(initBounce);
        if (n_ != n) {
            n_ = n;
            phi_.resize(n_);
            dr_ = rMax_ / (n_ - 1.);
            drinv_ = 1. / dr_;
            rToDimMin1_ = calcRToDimMin1();
        }
        initBounce.onGrid(phi_);
    }

    // return the number of the grid
    std::size_t n() const { return n_; }

    // return the dimension of space
    std::size_t dim() const { return dim_; }

    // return the lattice spacing
    double dr() const { return dr_; }

    //! The values of \f$ r^{dim - 1} \f$ on the grid.
    const std::vector<double> &rToDimMin1() const { return rToDimMin1_; }

    //! Compute the laplacian on the grid in radial coordinates. \f$\nabla^2\phi
    //! = \frac{d^2 phi}{dr^2}+(dim-1)/r*\frac{d phi}{dr} \f$.
    //!
    //! This function is time critical and therefore highly optimized. The
    //! laplacian for grid index zero is estimated from only two points, all
    //! other entries use three points (the last entry is never needed and thus
    //! not computed at all).
    std::vector<std::array<double, nPhi>> laplacian() const {
        std::vector<std::array<double, nPhi>> result(n_ - 1);
        const double drinvSq = std::pow(drinv_, 2);

        auto curr = phi_[0].begin();
        auto next = phi_[1].begin();
        auto out = result[0].begin();
        while (next != phi_[1].end()) {
            *out++ = 2 * (*next++ - *curr++) * drinvSq * dim_;
        }

        auto prev = phi_[0].begin();
        std::size_t elemCount = nPhi;
        const auto end = phi_.back().end();
        while (next != end) {
            const std::size_t i = elemCount++ / nPhi;
            const double dimFac = (dim_ - 1) / (2. * i);
            *out++ = drinvSq * ((1 + dimFac) * *next++ - 2 * *curr++ +
                                (1 - dimFac) * *prev++);
        }
        return result;
    }

    template <class Model>
    std::vector<std::array<double, nPhi>> dVdphi(const Model &model) {
        static_assert(Model::nPhi == nPhi, "field dimensions must match");
        std::vector<std::array<double, nPhi>> result(n_);
        auto phi = phi_.begin(), res = result.begin();
        const auto endRes = result.end();
        while (res != endRes) {
            model.calcDvdphi((phi++)->data(), (res++)->data());
        }
        return result;
    }

    // field excursion from the origin to the infinity
    double fieldExcursion() const {
        return Detail::normDifference(phi_.back(), phi_.front());
    }

    // derivative of scalar field at boundary
    double derivativeAtBoundary() const {
        return drinv_ * Detail::normDifference(phi_.back(), phi_[n_ - 2]);
    }
};

// kinetic energy of the configuration
// \int_0^\infty dr r^{d-1} \sum_i (-1/2) \phi_i \nabla^2\phi_i
template <std::size_t nPhi>
double kineticEnergy(const FieldConfiguration<nPhi> &field) {
    auto laplacian{field.laplacian()};
    auto integrand{field.rToDimMin1()};
    for (std::size_t i = 0; i != field.n() - 1; ++i) {
        integrand[i] *=
            -0.5 * std::inner_product(field[i].begin(), field[i].end(),
                                      laplacian[i].begin(), 0.);
    }
    return Detail::trapezoidalIntegrate(integrand.begin(), integrand.end(),
                                        field.dr());
}

// potential energy of the configuration
// \int_0^\infty dr r^{d-1} V(\phi)
template <class Model>
double potentialEnergy(const FieldConfiguration<Model::nPhi> &field,
                       const Model &model, double zeroPot) {
    auto integrand{field.rToDimMin1()};
    for (std::size_t i = 0; i != field.n(); ++i) {
        integrand[i] *= (model.vpot(field[i].data()) - zeroPot);
    }
    return Detail::trapezoidalIntegrate(integrand.begin(), integrand.end(),
                                        field.dr());
}

template <class Model> class BounceCalculator {
  public:
    static constexpr std::size_t nPhi = Model::nPhi;

  private:
    Parameters params_;
    std::size_t dim_;

    double lambda;

    const Model &model_;
    double VFV;

    //! Calculate \f$\lambda\f$ according to Eq. (5) of 1908.10868
    double computeLambda(const std::vector<std::array<double, nPhi>> &laplacian,
                         const std::vector<std::array<double, nPhi>> &dVdphi,
                         const std::vector<double> &rToDimMin1,
                         double gridSpacing) {
        auto integrand1{rToDimMin1};
        auto integrand2{rToDimMin1};

        auto i1{integrand1.begin()}, i2{integrand2.begin()};
        auto lap{laplacian.begin()}, dV{dVdphi.begin()};
        while (i1 != integrand1.end()) {
            *i1++ *=
                std::inner_product(dV->begin(), dV->end(), lap->begin(), 0.);
            *i2++ *=
                std::accumulate(dV->begin(), dV->end(), 0., Detail::squaredSum);
            ++lap;
            ++dV;
        }
        return Detail::trapezoidalIntegrate(integrand1.begin(),
                                            integrand1.end(), gridSpacing) /
               Detail::trapezoidalIntegrate(integrand2.begin(),
                                            integrand2.end(), gridSpacing);
    };

    //! Calculate the bracketed part of the RHS of Eq. (10) of 1908.10868.
    std::vector<std::array<double, nPhi>>
    computeRhs(std::vector<std::array<double, nPhi>> &&laplacian,
               const std::vector<std::array<double, nPhi>> &dVdphi,
               double lambda) {
        auto iRhs{laplacian.front().begin()};
        const auto endRhs{laplacian.back().end()};
        auto dV{dVdphi.front().begin()};
        while (iRhs != endRhs) {
            *iRhs++ -= lambda * *dV++;
        }
        return std::move(laplacian);
    }

    //! Perform a flow step for the field as in Eq. (10) of 1908.10868.
    //! The field value at boundary is fixed to phiFV and will not be updated.
    //! If the RHS of the EOM at the origin is too big a smaller step is taken.
    void performFlowStep(FieldConfiguration<nPhi> &field, double dTau,
                         const std::vector<std::array<double, nPhi>> &RHS) {
        const auto normAtOrigin{std::sqrt(std::accumulate(
            RHS.front().begin(), RHS.front().end(), 0., Detail::squaredSum))};
        const auto dtautilde{
            std::min(dTau, params_.maximumvariation * field.fieldExcursion() /
                               normAtOrigin)};

        auto iRhs{RHS.front().begin()};
        const auto endRhs{RHS.back().end()};
        auto f{field[0].begin()};
        while (iRhs != endRhs) {
            *f++ += dtautilde * *iRhs++;
        }
    }

    //! Euclidean area/volume term that appears in the bounce action.
    double bounceArea() {
        return dim_ * std::pow(M_PI, dim_ / 2.) / std::tgamma(dim_ / 2. + 1.);
    }

  public:
    BounceCalculator(const Model &model, std::size_t dim,
                     Parameters params = {})
        : params_{params}, dim_{dim}, model_{model} {}

    // main routine to get the bounce solution
    // See Fig. 1 of the manual
    double solve(const std::array<double, nPhi> &phiFV,
                 const std::array<double, nPhi> &phiTV) {
        VFV = model_.vpot(phiFV.data());
        if (model_.vpot(phiTV.data()) > VFV) {
            throw(std::invalid_argument(
                "The false vacuum is deeper than the true vacuum"));
        }

#ifdef SIMPLEBOUNCE_VERBOSE
        std::cerr << "probing a thickness to get negative V[phi] ..."
                  << std::endl;
#endif
        auto initBounce = InitialBounceConfiguration<nPhi>{
            phiFV, phiTV, params_.width, params_.r0Frac};
        // Make the bubble wall thin enough to get a negative potential energy.
        while (!initBounce.negativePotentialEnergy(model_, dim_)) {
            initBounce.width() *= 0.5;
        }

#ifdef SIMPLEBOUNCE_VERBOSE
        std::cerr << "probing the size of the bounce configuration ..."
                  << std::endl;
#endif
        FieldConfiguration<nPhi> field{initBounce, dim_, params_.rMax};
        // Make the size of the bubble sufficiently smaller than the size of the
        // sphere. If dphi/dr at the boundary becomes too large during flow,
        // take a smaller initial bounce configuration.
        while (evolveUntil(field,
                           params_.tend1 * params_.rMax * params_.rMax) != 0) {
#ifdef SIMPLEBOUNCE_VERBOSE
            std::cerr << "the size of the bounce is too large. initial "
                         "condition is scale transformed."
                      << std::endl;
#endif
            initBounce.r0byR() *= 0.5;
            initBounce.width() *= 0.5;
            field.resetInitialBounce(initBounce);
        }
        return action(field);
    }

    // evolve the configuration from tau = 0 to tau = tauend
    int evolveUntil(FieldConfiguration<nPhi> &field, const double tauend) {

        // 1 + d + sqrt(1 + d) is maximum of absolute value of eigenvalue of
        // {{-2d, 2d},{ (1-d)/2 + 1, -2}}, which is discreitzed Laplacian
        // for n = 2. This value is 6 for d=3, and 7.23607 for d=4. The
        // numerical value of maximum of absolute value of eigenvalue of
        // discretized Laplacian for large n is 6 for d=3, and 7.21417 for
        // d=4
        double dtau = 2. / (1. + dim_ + sqrt(1. + dim_)) * pow(field.dr(), 2) *
                      params_.safetyfactor;

#ifdef SIMPLEBOUNCE_VERBOSE
        std::cerr << "evolve until tau = " << tauend << ", (dtau = " << dtau
                  << ")" << std::endl;
#endif

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
        auto laplacian{field.laplacian()};
        const auto dVdphi{field.dVdphi(model_)};
        lambda =
            computeLambda(laplacian, dVdphi, field.rToDimMin1(), field.dr());
        const auto RHS = computeRhs(std::move(laplacian), dVdphi, lambda);
        performFlowStep(field, dtau, RHS);
        return lambda;
    }

    //! Kinetic energy of the bounce.
    double tBounce(const FieldConfiguration<nPhi> &field) const {
        return bounceArea() * std::pow(lambda, dim_ / 2. - 1.) *
               kineticEnergy(field);
    }

    //! Potential energy of the bounce.
    double vBounce(const FieldConfiguration<nPhi> &field) const {
        return bounceArea() * pow(lambda, dim_ / 2.) *
               potentialEnergy(field, model_, VFV);
    }

    //! Euclidean action of the bounce solution in d-dimensional space.
    double action(const FieldConfiguration<nPhi> &field) const {
        return tBounce(field) + vBounce(field);
    }
};

template <class Model, typename... Args>
auto makeBounceCalculator(const Model &model, std::size_t dim, Args... args) {
    return BounceCalculator<Model>{model, dim, std::forward<Args>(args)...};
}

} // namespace simplebounce
