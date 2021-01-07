#include "SimpleBounce/SimpleBounce.hpp"
#include "TestModels.hpp"
#include "catch.hpp"
namespace {

template <std::size_t nPhi> class OldFieldConfiguration {
    std::size_t n_, dim_;
    double rMax_, dr_, drinv_;
    double *phi_;
    std::vector<double> r_dminusoneth_;

    // radius : r_i = i * dr
    double r(int i) const { return dr_ * i; }

    // calculate and save pow(r(i), dim-1)
    void updateInfo() {
        r_dminusoneth_.resize(n_);
        for (int i = 0; i < n_; i++) {
            r_dminusoneth_[i] = pow(r(i), dim_ - 1);
        }
    }

    double &phi(std::size_t i, std::size_t iPhi) {
        return phi_[i * nPhi + iPhi];
    }
    double phi(std::size_t i, std::size_t iPhi) const {
        return phi_[i * nPhi + iPhi];
    }

  public:
    double *operator[](std::size_t i) { return &phi_[i * nPhi]; }
    const double *operator[](std::size_t i) const { return &phi_[i * nPhi]; }

    OldFieldConfiguration(std::size_t n, std::size_t dim, double rMax)
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
                phi(i, iphi) =
                    phiTV[iphi] + (phiFV[iphi] - phiTV[iphi]) *
                                      (1. + tanh((i - (n_ - 1) * frac) /
                                                 ((n_ - 1) * width))) /
                                      2.;
            }
        }
        for (int iphi = 0; iphi < nPhi; iphi++) {
            phi(n_ - 1, iphi) = phiFV[iphi];
        }
    }

    ~OldFieldConfiguration() { delete[] phi_; }

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
} // namespace

TEST_CASE("initialize field configuration") {
    static constexpr std::size_t nPhi = 5;
    const auto fv = std::array<double, nPhi>{0, 0, 0, 0, 0};
    const auto tv = std::array<double, nPhi>{1, 1, 1, 1, 1};
    BENCHMARK("old fieldconfig") {
        OldFieldConfiguration<nPhi> field(100, 4, 1.);
        field.setInitial(0.5, 0.05, fv, tv);
        return field;
    };

    BENCHMARK("new fieldconfig") {
        auto initField =
            simplebounce::InitialBounceConfiguration<nPhi>{fv, tv, 0.05, 0.5};
        auto field = simplebounce::FieldConfiguration<nPhi>{initField, 4, 1.};
        return field;
    };
}

TEST_CASE("initial configuration") {
    static constexpr std::size_t nPhi = 2;

    const auto fv = std::array<double, nPhi>{0, 0};
    const auto tv = std::array<double, nPhi>{1, 1};
    OldFieldConfiguration<nPhi> field{100, 4, 1.};
    field.setInitial(0.5, 0.05, fv, tv);
    simplebounce::InitialBounceConfiguration<nPhi> initField{fv, tv, 0.05, 0.5};

    auto getField = [&field](std::size_t i) {
        auto res = std::array<double, nPhi>{};
        std::copy(field[i], field[i] + nPhi, res.begin());
        return res;
    };
    CHECK(getField(0) == initField.phi(0));
    CHECK(getField(50)[0] == Approx(initField.phi(50 / 99.)[0]));
    CHECK(getField(60)[0] == Approx(initField.phi(60 / 99.)[0]));
    CHECK(getField(99) == initField.phi(1.));
    auto grid = initField.onGrid(100);
    CHECK(getField(0)[0] == Approx(grid[0][0]));
    CHECK(getField(99)[0] == Approx(grid[99][0]));
    CHECK(getField(40)[0] == Approx(grid[40][0]));
    CHECK(getField(50)[0] == Approx(grid[50][0]));

    auto field2 = simplebounce::FieldConfiguration<nPhi>{initField, 4, 1.};

    auto getField2 = [&field](std::size_t i) {
        auto res = std::array<double, nPhi>{};
        std::copy(field[i], field[i] + nPhi, res.begin());
        return res;
    };

    CHECK(getField2(0) == initField.phi(0));
    CHECK(getField2(50)[0] == Approx(initField.phi(50 / 99.)[0]));
    CHECK(getField2(60)[0] == Approx(initField.phi(60 / 99.)[0]));
    CHECK(getField2(99) == initField.phi(1.));

    CHECK(getField2(0) == getField(0));
    CHECK(getField2(50) == getField(50));
    CHECK(getField2(60) == getField(60));
    CHECK(getField2(99) == getField(99));
}

TEST_CASE("laplacian") {
    static constexpr std::size_t nPhi = 2;
    const auto fv = std::array<double, nPhi>{0, 0};
    const auto tv = std::array<double, nPhi>{1, 1};
    simplebounce::InitialBounceConfiguration<nPhi> initField{fv, tv, 0.05, 0.5};
    simplebounce::FieldConfiguration<nPhi> field{initField, 4, 1.};

    BENCHMARK("laplacian new") { return field.laplacian(); };

    OldFieldConfiguration<nPhi> oldField{100, 4, 1.};
    oldField.setInitial(0.5, 0.05, fv, tv);

    BENCHMARK("laplacian old") {
        const auto n = oldField.n();
        std::vector<std::array<double, nPhi>> laplacian(n);

        // \nabla^2 \phi_iphi at r = r_i
        for (int i = 0; i < n - 1; i++) {
            for (int iphi = 0; iphi < nPhi; iphi++) {
                laplacian[i][iphi] = oldField.lap(i, iphi);
            }
        }
        return laplacian;
    };
}

TEST_CASE("evolution") {

    static constexpr std::size_t nPhi = 2;
    static constexpr std::size_t dim = 4ul;
    const auto fv = std::array<double, nPhi>{0, 0};
    const auto tv = std::array<double, nPhi>{1, 1};
    simplebounce::InitialBounceConfiguration<nPhi> initField{fv, tv, 0.05, 0.5};
    simplebounce::FieldConfiguration<nPhi> field{initField, 4, 1.};

    double dtau = 2. / (1. + dim + sqrt(1. + dim)) * pow(field.dr(), 2) * 0.9;

    auto model = Model2();
    BENCHMARK("first step") {
        const auto n = field.n();
        auto laplacian = field.laplacian();

        // integral1 : \int_0^\infty dr r^{d-1} \sum_i (\partial V /
        // \partial\phi_i) \nabla^2\phi_i integral2 : \int_0^\infty dr
        // r^{d-1} \sum_i (\partial V / \partial\phi_i)^2
        auto integrand1{field.rToDimMin1()};
        auto integrand2{field.rToDimMin1()};
        for (int i = 0; i != n - 1; ++i) {
            std::array<double, nPhi> dvdphi;
            model.calcDvdphi(field[i], dvdphi.data());
            auto int1Sum{0.};
            for (int iphi = 0; iphi < nPhi; iphi++) {
                int1Sum += dvdphi[iphi] * laplacian[i][iphi];
            }
            integrand1[i] *= int1Sum;
            integrand2[i] *= std::accumulate(
                dvdphi.begin(), dvdphi.end(), 0.,
                [](double sum, double x) { return sum + std::pow(x, 2); });
        }

        // Eq. 9 of 1907.02417
        auto lambda = simplebounce::trapezoidalIntegrate(
                          integrand1.begin(), integrand1.end(), field.dr()) /
                      simplebounce::trapezoidalIntegrate(
                          integrand2.begin(), integrand2.end(), field.dr());

        // RHS of Eq. 8 of 1907.02417
        // phi at boundary is fixed to phiFV and will not be updated.
        auto RHS = std::move(laplacian);
        for (int i = 0; i < n - 1; i++) {
            std::array<double, nPhi> dvdphi;
            model.calcDvdphi(field[i], dvdphi.data());
            for (int iphi = 0; iphi < nPhi; iphi++) {
                RHS[i][iphi] -= lambda * dvdphi[iphi];
            }
        }

        // if RHS of EOM at the origin is too big, smaller step is taken.
        double sum = 0.;
        for (int iphi = 0; iphi < nPhi; iphi++) {
            sum += RHS[0][iphi] * RHS[0][iphi];
        }
        double dtautilde = 1e-2 * field.fieldExcursion() / sqrt(sum);
        if (dtau < dtautilde) {
            dtautilde = dtau;
        }
        return lambda;
    };
}
