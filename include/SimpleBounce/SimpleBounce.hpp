/** @file SimpleBounce.hpp
 *
 * Simplebounce is a single-header C++14 library that computes the bounce
 * configuration of vacuum tunnelling for models with arbitrary scalar
 * potentials. It is based on the method proposed in arXiv:1907.02417 and the
 * algorithm is documented in arXiv:1908.10868. All equation numbers mentioned
 * in this file refer to the latter publication.
 */
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

//! Namespace of the simplebounce library
namespace simplebounce {
template <class Model> struct BounceSolution;
template <std::size_t nPhi> class FieldConfiguration;
namespace Detail {
template <std::size_t nPhi> class InitialBounceConfiguration;
}

//! Parameters that control the simplebounce algorithm.
struct Parameters {
    // minimal grid size
    std::size_t minGridSize = 100;
    // radius \f$R\f$ of the boundary
    double rMax = 1.;
    // safety factor to ensure no overly large $\delta\tau$ step is taken
    double safetyfactor = 0.9;
    // maximal variation of the field in one evolution step relative to the
    // total field excursion
    double maximumvariation = 0.01;
    // position \f$r_0/R\f$ for the steep region of the initial bounce Eq (11)
    double r0byR = 0.5;
    // width \f$\sigma/R\f$ for the initial bounce Eq (11)
    double width = 0.05;
    // maximal value of the derivative of the field at the boundary
    double derivMax = 1e-2;
    // flow time \f$\tau_1\f$
    double tau1 = 0.4;
};

/**
 * Main routine of simplebounce. Calculates the bounce solution from the false
 * vacuum to the true vacuum.
 *
 * @tparam Model a Model class for use with simplebounce has to provide the
 * following member variables and functions:
 * ```c++
 * class ExampleModel {
 *     // number of field degrees of freedom
 *     static constexpr std::size_t nPhi = 2;
 *     // value of the scalar potential at the given fieldspace point phi
 *     // (with nPhi elements)
 *     double vpot(const double phi[]) const;
 *     // compute the derivative of the scalar potential at the fieldspace
 *     // point phi (both with nPhi elements)
 *     void calcDvdphi(const double phi[], double derivative[]);
 * };
 * ```
 * @param model Instance of a Model class that will be used to calculate the
 * bounce
 * @param phiFV Fieldspace location of the false vacuum
 * @param phiTV Fieldspace location of the true vacuum. Any point in the
 * potential deeper than the false vacuum is sufficient, this need not be a
 * minimum.
 * @param dim Dimensionality of the euclidean space in which the tunnelling
 * takes place. This should usually be either 4 for zero temperature, or 3 for
 * finite temperature tunnelling processes.
 * @param params Optionally adjust the parameters of the algorithm to your
 * liking.
 * @returns BounceSolution<Model> the bounce solution from phiFV to phiTV
 */
template <class Model>
BounceSolution<Model> solve(const Model &model,
                            const std::array<double, Model::nPhi> &phiFV,
                            const std::array<double, Model::nPhi> &phiTV,
                            std::size_t dim, Parameters params = {}) {
    if (model.vpot(phiTV.data()) > model.vpot(phiFV.data())) {
        throw(std::invalid_argument(
            "The false vacuum is deeper than the true vacuum"));
    }

#ifdef SIMPLEBOUNCE_VERBOSE
    std::cerr << "probing a thickness to get negative V[phi] ..." << std::endl;
#endif
    auto initBounce{Detail::InitialBounceConfiguration<Model::nPhi>{
        phiFV, phiTV, params.width, params.r0byR}};
    // Make the bubble wall thin enough to get a negative potential energy.
    while (!initBounce.negativePotentialEnergy(model, dim)) {
        initBounce.width() *= 0.5;
    }

#ifdef SIMPLEBOUNCE_VERBOSE
    std::cerr << "probing the size of the bounce configuration ..."
              << std::endl;
#endif
    BounceSolution<Model> result{
        {initBounce, dim, params.rMax, params.minGridSize}, model};
    // Make the size of the bubble sufficiently smaller than the size of the
    // sphere. If dphi/dr at the boundary becomes too large during flow,
    // take a smaller initial bounce configuration.
    while (!evolveUntil(result, params)) {
#ifdef SIMPLEBOUNCE_VERBOSE
        std::cerr << "the size of the bounce is too large. initial "
                     "condition is scale transformed."
                  << std::endl;
#endif
        initBounce.r0byR() *= 0.5;
        initBounce.width() *= 0.5;
        result.field.resetInitialBounce(initBounce);
    }
    return result;
}

//! A bounce solution calculated by simplebounce.
template <class Model> struct BounceSolution {
    //! The field values of the reduced bounce solution.
    FieldConfiguration<Model::nPhi> field;
    //! The physics model in which this solution was calculated.
    const Model &model;
    //! Scale transformation parameter \f$\lambda\f$ to obtain the bounce
    //! solution from the reduces bounce solution through Eq 8.
    double lambda = 0.;

    //! Euclidean action of the bounce solution in d-dimensional space.
    double action() const
        noexcept(noexcept(model.vpot(std::declval<double *>()))) {
        using std::pow;
        const auto do2{field.dim() / 2.};
        const auto area{2 * do2 * pow(M_PI, do2) / std::tgamma(do2 + 1.)};
        const auto tBounce{area * pow(lambda, do2 - 1.) * kineticEnergy(field)};
        const auto vBounce{area * pow(lambda, do2) *
                           potentialEnergy(field, model)};
        return tBounce + vBounce;
    }
};

//! Simplebounce implementation details.
namespace Detail {
//! Functor that performs a squared sum. E.g. for use with std::accumulate.
const auto squaredSum = [](double sum, double x) noexcept {
    return sum + std::pow(x, 2);
};

//! Computes the euclidean vector norm \f$ | \vec{x1} - \vec{x2} | \f$.
template <std::size_t n>
double normDifference(const std::array<double, n> &x1,
                      const std::array<double, n> &x2) noexcept {
    return std::sqrt(std::inner_product(x1.begin(), x1.end(), x2.begin(), 0.,
                                        Detail::squaredSum,
                                        std::minus<double>{}));
}

//! 1D trapezoidal integration on a uniform grid.
template <class BidirectionalIterator>
double trapezoidalIntegrate(BidirectionalIterator beginValues,
                            BidirectionalIterator endValues,
                            double gridSpacing) noexcept {
    double borderVal = (*beginValues++ + *endValues--) / 2.;
    return gridSpacing * std::accumulate(beginValues, endValues, borderVal);
}
} // namespace Detail

//! A field configuration in nPhi field dimensions on a uniformly spaced grid
//! in \f$r\f$ space.
template <std::size_t nPhi> class FieldConfiguration {
  public:
    //! maximal grid size
    static constexpr std::size_t maxN = 1000;

  private:
    std::size_t dim_;
    double rMax_;
    std::vector<std::array<double, nPhi>> phi_;

    double dr_ = rMax_ / (phi_.size() - 1);
    std::vector<double> rToDimMin1_ = calcRToDimMin1();

    // calculate \f$ r_i^{dim-1} \f$
    std::vector<double> calcRToDimMin1() const noexcept {
        std::vector<double> result(phi_.size(), std::pow(dr_, dim_ - 1));
        std::size_t i{0};
        for (auto &r : result) {
            r *= std::pow(i++, dim_ - 1);
        }
        return result;
    }

    std::size_t determineGridSize(
        const Detail::InitialBounceConfiguration<nPhi> &initBounce) const {
        const auto n{std::max(std::size_t{100}, static_cast<std::size_t>(
                                                    1 / initBounce.width()))};
        if (n > maxN) {
            throw(std::runtime_error("Maximum grid size exceeded"));
        }
        return n;
    }

  public:
    //! access the field values at grid index i
    std::array<double, nPhi> &operator[](std::size_t i) noexcept {
        return phi_[i];
    }
    //! access the field values at grid index i
    const std::array<double, nPhi> &operator[](std::size_t i) const noexcept {
        return phi_[i];
    }
    //! field values of the false vacuum
    const std::array<double, nPhi> &phiFV() const noexcept {
        return phi_.back();
    }
    //! begin iterator through the grid points
    auto begin() const noexcept { return phi_.begin(); }
    //! end iterator through the grid points
    auto end() const noexcept { return phi_.end(); }

    //! Construct a field configuration from an initial bounce configuration.
    FieldConfiguration(
        const Detail::InitialBounceConfiguration<nPhi> &initBounce,
        std::size_t dim, double rMax, std::size_t minGridSize)
        : dim_{dim}, rMax_{rMax}, phi_{initBounce.onGrid(
                                      determineGridSize(initBounce))} {}

    //! Reset to a different initial bounce configuration.
    void resetInitialBounce(
        const Detail::InitialBounceConfiguration<nPhi> &initBounce) {
        auto n{determineGridSize(initBounce)};
        if (phi_.size() != n) {
            phi_.resize(n);
            dr_ = rMax_ / (n - 1.);
            rToDimMin1_ = calcRToDimMin1();
        }
        initBounce.onGrid(phi_);
    }

    //! The dimension of euclidan space associated with this field configuration
    std::size_t dim() const noexcept { return dim_; }

    //! The grid spacing
    double dr() const noexcept { return dr_; }

    //! The values of \f$ r^{dim - 1} \f$ on the grid. This function exists only
    //! to speed up computation since these values are frequently needed.
    const std::vector<double> &rToDimMin1() const noexcept {
        return rToDimMin1_;
    }

    //! Compute the laplacian on the grid in radial coordinates. \f$\nabla^2\phi
    //! = \frac{d^2 phi}{dr^2}+(dim-1)/r*\frac{d phi}{dr} \f$.
    //!
    //! This function is time critical and therefore highly optimized. The
    //! laplacian for grid index zero is estimated from only two points, all
    //! other entries use three points (the last entry is never needed and thus
    //! not computed at all).
    std::vector<std::array<double, nPhi>> laplacian() const noexcept {
        std::vector<std::array<double, nPhi>> result(phi_.size() - 1);
        const double drinvSq{std::pow(dr_, -2)};

        auto curr{phi_[0].begin()};
        auto next{phi_[1].begin()};
        auto out{result[0].begin()};
        while (next != phi_[1].end()) {
            *out++ = 2 * (*next++ - *curr++) * drinvSq * dim_;
        }

        auto prev{phi_[0].begin()};
        std::size_t elemCount{nPhi};
        const auto end{phi_.back().end()};
        while (next != end) {
            const std::size_t i{elemCount++ / nPhi};
            const double dimFac{(dim_ - 1) / (2. * i)};
            *out++ = drinvSq * ((1 + dimFac) * *next++ - 2 * *curr++ +
                                (1 - dimFac) * *prev++);
        }
        return result;
    }

    //! Compute the gradient \f$d V(\vec\phi)/d\phi_i\f$ of the scalar potential
    //! of the model at every grid point.
    template <class Model>
    std::vector<std::array<double, nPhi>> dVdphi(const Model &model) const
        noexcept(noexcept(model.calcDvdphi(std::declval<double *>(),
                                           std::declval<double *>()))) {
        static_assert(Model::nPhi == nPhi, "field dimensions must match");
        std::vector<std::array<double, nPhi>> result(phi_.size());
        auto phi{phi_.begin()};
        for (auto &res : result) {
            model.calcDvdphi((phi++)->data(), res.data());
        }
        return result;
    }

    //! Total field excursion from \f$r=0\f$ to \f$r=R\f$.
    double fieldExcursion() const noexcept {
        return Detail::normDifference(phi_.back(), phi_.front());
    }

    //! Norm of the derivative of scalar field at boundary.
    double derivativeAtBoundary() const noexcept {
        return Detail::normDifference(phi_.back(), *(phi_.end() - 2)) / dr_;
    }
};

//! Kinetic energy of the field configuration.
//! \f$ \int_0^\infty dr r^{d-1} (-1/2) \sum_i \phi_i \nabla^2\phi_i \f$
template <std::size_t nPhi>
double kineticEnergy(const FieldConfiguration<nPhi> &field) noexcept {
    auto laplacian{field.laplacian()};
    auto integrand{field.rToDimMin1()};
    auto iField{field.begin()};
    auto iInteg{integrand.begin()};
    for (const auto &lap : laplacian) {
        *iInteg++ *= -0.5 * std::inner_product(lap.begin(), lap.end(),
                                               (iField++)->begin(), 0.);
    }
    return Detail::trapezoidalIntegrate(integrand.begin(), integrand.end(),
                                        field.dr());
}

//! Potential energy of the field configuration in the given model.
//! \f$ \int_0^\infty dr r^{d-1} V(\phi) \f$
template <class Model>
double potentialEnergy(
    const FieldConfiguration<Model::nPhi> &field,
    const Model
        &model) noexcept(noexcept(model.vpot(std::declval<double *>()))) {
    auto integrand{field.rToDimMin1()};
    auto iField{field.begin()};
    auto zeroPot{model.vpot(field.phiFV().data())};
    for (double &integ : integrand) {
        integ *= (model.vpot((iField++)->data()) - zeroPot);
    }
    return Detail::trapezoidalIntegrate(integrand.begin(), integrand.end(),
                                        field.dr());
}

namespace Detail {

//! Implements the initial bounce configuration as defined in Eq (11). This
//! class allows adjusting the initial configuration without having to construct
//! the whole grid of a FieldConfiguration.
template <std::size_t nPhi> class InitialBounceConfiguration {
    std::array<double, nPhi> phiFV_;
    std::array<double, nPhi> phiTV_;
    double width_;
    double r0byR_;

    static constexpr double half = 1 / 2.;
    static constexpr double almost0 = 0.1 / FieldConfiguration<nPhi>::maxN;
    static constexpr double almost1 = 1 - almost0;

  public:
    //! Construct an initial configuration between phiFV and phiTV according to
    //! Eq (11).
    InitialBounceConfiguration(const std::array<double, nPhi> &phiFV,
                               const std::array<double, nPhi> &phiTV,
                               double width, double r0byR) noexcept
        : phiFV_{phiFV}, phiTV_{phiTV}, width_{width}, r0byR_{r0byR} {}

    //! Value of the field configuration at the given value of \f$r/R\f$. By
    //! construction, the false vacuum is reached as rbyR->1 and the true vacuum
    //! as rbyR->0.
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

    //! Construct a grid of \f$n\f$ evenly spaced points of phi in rbyR.
    std::vector<std::array<double, nPhi>> onGrid(std::size_t n) const noexcept {
        auto phiGrid{std::vector<std::array<double, nPhi>>(n)};
        onGrid(phiGrid);
        return phiGrid;
    }

    //! Fill the given grid with evenly spaced points of phi in rbyR.
    void onGrid(std::vector<std::array<double, nPhi>> &grid) const noexcept {
        const auto n{grid.size()};
        grid.front() = phiTV_;
        const auto delta{1. / static_cast<double>(n - 1)};
        for (std::size_t i = 1; i != n - 1; ++i) {
            grid[i] = phi(i * delta);
        }
        grid.back() = phiFV_;
    }

    double width() const noexcept { return width_; } //!< \f$\sigma/R\f$
    double &width() noexcept { return width_; }      //!< \f$\sigma/R\f$
    double &r0byR() noexcept { return r0byR_; }      //!< \f$r_0/R\f$

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
    bool negativePotentialEnergy(const Model &model, std::size_t dim) const
        noexcept(noexcept(model.vpot(std::declval<double *>()))) {
        static constexpr std::size_t maxRefinements = 10;
        static constexpr double margin = 5;
        static constexpr double maxVal = -1e-5;

        const auto zeroPot{model.vpot(phiFV_.data())};
        const auto integrand = [dim, zeroPot, model, this](double r) {
            return std::pow(r, dim - 1) * (model.vpot(phi(r).data()) - zeroPot);
        };

        // the integrand vanishes at the endpoints by construction
        auto I0{0.};
        auto h{half};
        auto I1{integrand(h) * h};

        // The recursion is:
        // I_k = 1/2 I_{k-1} + 1/2^k \sum_{j=1; j odd, j < 2^k} f(a+j(b-a)/2^k)
        auto k{2ul};
        auto error{std::abs(I0 - I1)};
        while (k < 4 || (k < maxRefinements && margin * error > std::abs(I1))) {
            I0 = I1;
            I1 = I0 * half;
            const auto p{static_cast<std::size_t>(1) << k};
            h *= half;
            auto sum{0.};
            for (std::size_t j = 1; j < p; j += 2) {
                auto y{integrand(j * h)};
                sum += y;
            }
            I1 += sum * h;
            ++k;
            error = std::abs(I0 - I1);
        }
        return I1 < maxVal - error;
    }
};

//! Calculate \f$\lambda\f$ according to Eq (5).
template <std::size_t nPhi>
double computeLambda(const std::vector<std::array<double, nPhi>> &laplacian,
                     const std::vector<std::array<double, nPhi>> &dVdphi,
                     const std::vector<double> &rToDimMin1,
                     double gridSpacing) noexcept {
    auto integrand1{rToDimMin1};
    auto integrand2{rToDimMin1};

    auto i1{integrand1.begin()}, i2{integrand2.begin()};
    auto dV{dVdphi.begin()};
    for (const auto &lap : laplacian) {
        *i1++ *= std::inner_product(dV->begin(), dV->end(), lap.begin(), 0.);
        *i2++ *=
            std::accumulate(dV->begin(), dV->end(), 0., Detail::squaredSum);
        ++dV;
    }
    integrand1.back() = 0;
    integrand2.back() = 0;
    return Detail::trapezoidalIntegrate(integrand1.begin(), integrand1.end(),
                                        gridSpacing) /
           Detail::trapezoidalIntegrate(integrand2.begin(), integrand2.end(),
                                        gridSpacing);
};

//! Calculate the bracketed part of the RHS of Eq (10).
//! Reuse the memory from laplacian as an optimization.
template <std::size_t nPhi>
std::vector<std::array<double, nPhi>>
computeRhs(std::vector<std::array<double, nPhi>> &&laplacian,
           const std::vector<std::array<double, nPhi>> &dVdphi,
           double lambda) noexcept {
    auto iRhs{laplacian.front().begin()};
    const auto endRhs{laplacian.back().end()};
    auto dV{dVdphi.front().begin()};
    while (iRhs != endRhs) {
        *iRhs++ -= lambda * *dV++;
    }
    return std::move(laplacian);
}

//! Perform a flow step for the field as in Eq (10) of 1908.10868.
//! The field value at upper boundary is fixed to phiFV and will not be updated.
//! If the RHS of the EOM at the origin is too big a smaller step is taken.
template <std::size_t nPhi>
void performFlowStep(FieldConfiguration<nPhi> &field, double dTau,
                     const std::vector<std::array<double, nPhi>> &RHS,
                     double maxVariation) noexcept {
    const auto normAtOrigin{std::sqrt(std::accumulate(
        RHS.front().begin(), RHS.front().end(), 0., Detail::squaredSum))};
    const auto dtautilde{
        std::min(dTau, maxVariation * field.fieldExcursion() / normAtOrigin)};

    auto iRhs{RHS.front().begin()};
    const auto endRhs{RHS.back().end()};
    auto f{field[0].begin()};
    while (iRhs != endRhs) {
        *f++ += dtautilde * *iRhs++;
    }
}
} // namespace Detail

//! Evolve the field configuration by dtau as in Eq (10).
template <class Model>
void evolve(
    BounceSolution<Model> &solution, double dtau,
    double maxVariation) noexcept(noexcept(solution.field
                                               .dVdphi(solution.model))) {
    auto laplacian{solution.field.laplacian()};
    const auto dVdphi{solution.field.dVdphi(solution.model)};
    solution.lambda = Detail::computeLambda(
        laplacian, dVdphi, solution.field.rToDimMin1(), solution.field.dr());
    const auto RHS =
        Detail::computeRhs(std::move(laplacian), dVdphi, solution.lambda);
    Detail::performFlowStep(solution.field, dtau, RHS, maxVariation);
}

//! Evolve the field configuration from `tau = 0` to `tau = tauEnd`. Aborts and
//! returns false if boundary derivative becomes too large.
template <class Model>
bool evolveUntil(
    BounceSolution<Model> &solution,
    const Parameters
        &params) noexcept(noexcept(std::declval<Model>()
                                       .calcDvdphi(std::declval<double *>(),
                                                   std::declval<double *>()))) {

    const double tauEnd = params.tau1 * params.rMax * params.rMax;
    // 1 + d + sqrt(1 + d) is maximum of absolute value of eigenvalue of
    // {{-2d, 2d},{ (1-d)/2 + 1, -2}}, which is discreitzed Laplacian
    // for n = 2. This value is 6 for d=3, and 7.23607 for d=4. The
    // numerical value of maximum of absolute value of eigenvalue of
    // discretized Laplacian for large n is 6 for d=3, and 7.21417 for
    // d=4
    const std::size_t dim = solution.field.dim();
    const double dtau = 2. / (1. + dim + std::sqrt(1. + dim)) *
                        std::pow(solution.field.dr(), 2) * params.safetyfactor;

#ifdef SIMPLEBOUNCE_VERBOSE
    std::cerr << "evolve until tau = " << tauEnd << ", (dtau = " << dtau << ")"
              << std::endl;
#endif

    for (double tau = 0.; tau < tauEnd; tau += dtau) {
        evolve(solution, dtau, params.maximumvariation);
        if (solution.field.derivativeAtBoundary() * params.rMax /
                solution.field.fieldExcursion() >
            params.derivMax) {
            return false;
        }
    }
    return true;
}

} // namespace simplebounce
