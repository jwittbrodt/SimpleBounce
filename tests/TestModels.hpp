#include "SimpleBounce/SimpleBounce.hpp"

class Model1 : public simplebounce::GenericModel<1> {
  public:
    double vpot(const double *phi) const {
        return (phi[0] * phi[0] * phi[0] * phi[0] -
                8. * phi[0] * phi[0] * phi[0] + 10. * phi[0] * phi[0]) /
               10.;
    }
    void calcDvdphi(const double *phi, double *dvdphi) const {
        dvdphi[0] = 0.4 * phi[0] * phi[0] * phi[0] - 2.4 * phi[0] * phi[0] +
                    2. * phi[0];
    }
};

class Model1a : public simplebounce::GenericModel<1> {
  public:
    static constexpr double c = 0.47;
    double vpot(const double *phi) const {
        return phi[0] * phi[0] * phi[0] * phi[0] / 4. -
               (c + 1.) / 3. * phi[0] * phi[0] * phi[0] +
               c / 2. * phi[0] * phi[0];
    }
    void calcDvdphi(const double *phi, double *dvdphi) const {
        dvdphi[0] =
            phi[0] * phi[0] * phi[0] - (c + 1.) * phi[0] * phi[0] + c * phi[0];
    }
};

class Model1b : public simplebounce::GenericModel<1> {
  public:
    static constexpr double c = 0.2;
    double vpot(const double *phi) const {
        return phi[0] * phi[0] * phi[0] * phi[0] / 4. -
               (c + 1.) / 3. * phi[0] * phi[0] * phi[0] +
               c / 2. * phi[0] * phi[0];
    }
    void calcDvdphi(const double *phi, double *dvdphi) const {
        dvdphi[0] =
            phi[0] * phi[0] * phi[0] - (c + 1.) * phi[0] * phi[0] + c * phi[0];
    }
};

class Model2 : public simplebounce::GenericModel<2> {
  public:
    static constexpr double c0 = 1.8;
    static constexpr double c1 = 0.2;
    static constexpr double c2 = 0.3;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1]);
        return (r1 - c2) * r2;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1]);
        dvdphi[0] = 2. * c0 * (phi[0] - 1.) * r2 + 2. * phi[0] * (r1 - c2);
        dvdphi[1] = 2. * c1 * (phi[1] - 1.) * r2 + 2. * phi[1] * (r1 - c2);
    }
};

class Model2a : public simplebounce::GenericModel<2> {
  public:
    static constexpr double c = 2.;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = phi[0] * phi[0] + 5. * phi[1] * phi[1];
        double r2 =
            5. * (phi[0] - 1.) * (phi[0] - 1.) + (phi[1] - 1.) * (phi[1] - 1.);
        double r3 = c * (phi[1] * phi[1] * phi[1] * phi[1] / 4. -
                         phi[1] * phi[1] * phi[1] / 3.);
        return r1 * r2 + r3;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = phi[0] * phi[0] + 5. * phi[1] * phi[1];
        double r2 =
            5. * (phi[0] - 1.) * (phi[0] - 1.) + (phi[1] - 1.) * (phi[1] - 1.);
        double dr1dx = 2. * phi[0];
        double dr1dy = 10. * phi[1];
        double dr2dx = 10. * (phi[0] - 1.);
        double dr2dy = 2. * (phi[1] - 1.);
        double dr3dy = c * (phi[1] * phi[1] * phi[1] - phi[1] * phi[1]);
        dvdphi[0] = dr1dx * r2 + r1 * dr2dx;
        dvdphi[1] = dr1dy * r2 + r1 * dr2dy + dr3dy;
    }
};

class Model2b : public simplebounce::GenericModel<2> {
  public:
    static constexpr double c = 80.;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = phi[0] * phi[0] + 5. * phi[1] * phi[1];
        double r2 =
            5. * (phi[0] - 1.) * (phi[0] - 1.) + (phi[1] - 1.) * (phi[1] - 1.);
        double r3 = c * (phi[1] * phi[1] * phi[1] * phi[1] / 4. -
                         phi[1] * phi[1] * phi[1] / 3.);
        return r1 * r2 + r3;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = phi[0] * phi[0] + 5. * phi[1] * phi[1];
        double r2 =
            5. * (phi[0] - 1.) * (phi[0] - 1.) + (phi[1] - 1.) * (phi[1] - 1.);
        double dr1dx = 2. * phi[0];
        double dr1dy = 10. * phi[1];
        double dr2dx = 10. * (phi[0] - 1.);
        double dr2dy = 2. * (phi[1] - 1.);
        double dr3dy = c * (phi[1] * phi[1] * phi[1] - phi[1] * phi[1]);
        dvdphi[0] = dr1dx * r2 + r1 * dr2dx;
        dvdphi[1] = dr1dy * r2 + r1 * dr2dy + dr3dy;
    }
};

class Model3 : public simplebounce::GenericModel<3> {
  public:
    static constexpr double c0 = 0.684373;
    static constexpr double c1 = 0.181928;
    static constexpr double c2 = 0.295089;
    static constexpr double c3 = 0.284821;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);
        return (r1 - c3) * r2;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);
        dvdphi[0] = 2. * c0 * (phi[0] - 1.) * r2 + 2. * phi[0] * (r1 - c3);
        dvdphi[1] = 2. * c1 * (phi[1] - 1.) * r2 + 2. * phi[1] * (r1 - c3);
        dvdphi[2] = 2. * c2 * (phi[2] - 1.) * r2 + 2. * phi[2] * (r1 - c3);
    }
};

class Model4 : public simplebounce::GenericModel<4> {
  public:
    static constexpr double c0 = 0.534808;
    static constexpr double c1 = 0.77023;
    static constexpr double c2 = 0.838912;
    static constexpr double c3 = 0.00517238;
    static constexpr double c4 = 0.258889;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3]);
        return (r1 - c4) * r2;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3]);
        dvdphi[0] = 2. * c0 * (phi[0] - 1.) * r2 + 2. * phi[0] * (r1 - c4);
        dvdphi[1] = 2. * c1 * (phi[1] - 1.) * r2 + 2. * phi[1] * (r1 - c4);
        dvdphi[2] = 2. * c2 * (phi[2] - 1.) * r2 + 2. * phi[2] * (r1 - c4);
        dvdphi[3] = 2. * c3 * (phi[3] - 1.) * r2 + 2. * phi[3] * (r1 - c4);
    }
};

class Model5 : public simplebounce::GenericModel<5> {
  public:
    static constexpr double c0 = 0.4747;
    static constexpr double c1 = 0.234808;
    static constexpr double c2 = 0.57023;
    static constexpr double c3 = 0.138912;
    static constexpr double c4 = 0.517238;
    static constexpr double c5 = 0.658889;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.) +
                     c4 * (phi[4] - 1.) * (phi[4] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3] + phi[4] * phi[4]);
        return (r1 - c5) * r2;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.) +
                     c4 * (phi[4] - 1.) * (phi[4] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3] + phi[4] * phi[4]);
        dvdphi[0] = 2. * c0 * (phi[0] - 1.) * r2 + 2. * phi[0] * (r1 - c5);
        dvdphi[1] = 2. * c1 * (phi[1] - 1.) * r2 + 2. * phi[1] * (r1 - c5);
        dvdphi[2] = 2. * c2 * (phi[2] - 1.) * r2 + 2. * phi[2] * (r1 - c5);
        dvdphi[3] = 2. * c3 * (phi[3] - 1.) * r2 + 2. * phi[3] * (r1 - c5);
        dvdphi[4] = 2. * c4 * (phi[4] - 1.) * r2 + 2. * phi[4] * (r1 - c5);
    }
};

class Model6 : public simplebounce::GenericModel<6> {
  public:
    static constexpr double c0 = 0.34234;
    static constexpr double c1 = 0.4747;
    static constexpr double c2 = 0.234808;
    static constexpr double c3 = 0.57023;
    static constexpr double c4 = 0.138912;
    static constexpr double c5 = 0.517238;
    static constexpr double c6 = 0.658889;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.) +
                     c4 * (phi[4] - 1.) * (phi[4] - 1.) +
                     c5 * (phi[5] - 1.) * (phi[5] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3] + phi[4] * phi[4] + phi[5] * phi[5]);
        return (r1 - c6) * r2;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.) +
                     c4 * (phi[4] - 1.) * (phi[4] - 1.) +
                     c5 * (phi[5] - 1.) * (phi[5] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3] + phi[4] * phi[4] + phi[5] * phi[5]);
        dvdphi[0] = 2. * c0 * (phi[0] - 1.) * r2 + 2. * phi[0] * (r1 - c6);
        dvdphi[1] = 2. * c1 * (phi[1] - 1.) * r2 + 2. * phi[1] * (r1 - c6);
        dvdphi[2] = 2. * c2 * (phi[2] - 1.) * r2 + 2. * phi[2] * (r1 - c6);
        dvdphi[3] = 2. * c3 * (phi[3] - 1.) * r2 + 2. * phi[3] * (r1 - c6);
        dvdphi[4] = 2. * c4 * (phi[4] - 1.) * r2 + 2. * phi[4] * (r1 - c6);
        dvdphi[5] = 2. * c5 * (phi[5] - 1.) * r2 + 2. * phi[5] * (r1 - c6);
    }
};

class Model7 : public simplebounce::GenericModel<7> {
  public:
    static constexpr double c0 = 0.5233;
    static constexpr double c1 = 0.34234;
    static constexpr double c2 = 0.4747;
    static constexpr double c3 = 0.234808;
    static constexpr double c4 = 0.57023;
    static constexpr double c5 = 0.138912;
    static constexpr double c6 = 0.517238;
    static constexpr double c7 = 0.65889;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.) +
                     c4 * (phi[4] - 1.) * (phi[4] - 1.) +
                     c5 * (phi[5] - 1.) * (phi[5] - 1.) +
                     c6 * (phi[6] - 1.) * (phi[6] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3] + phi[4] * phi[4] + phi[5] * phi[5] +
                     phi[6] * phi[6]);
        return (r1 - c7) * r2;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.) +
                     c4 * (phi[4] - 1.) * (phi[4] - 1.) +
                     c5 * (phi[5] - 1.) * (phi[5] - 1.) +
                     c6 * (phi[6] - 1.) * (phi[6] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3] + phi[4] * phi[4] + phi[5] * phi[5] +
                     phi[6] * phi[6]);
        dvdphi[0] = 2. * c0 * (phi[0] - 1.) * r2 + 2. * phi[0] * (r1 - c7);
        dvdphi[1] = 2. * c1 * (phi[1] - 1.) * r2 + 2. * phi[1] * (r1 - c7);
        dvdphi[2] = 2. * c2 * (phi[2] - 1.) * r2 + 2. * phi[2] * (r1 - c7);
        dvdphi[3] = 2. * c3 * (phi[3] - 1.) * r2 + 2. * phi[3] * (r1 - c7);
        dvdphi[4] = 2. * c4 * (phi[4] - 1.) * r2 + 2. * phi[4] * (r1 - c7);
        dvdphi[5] = 2. * c5 * (phi[5] - 1.) * r2 + 2. * phi[5] * (r1 - c7);
        dvdphi[6] = 2. * c6 * (phi[6] - 1.) * r2 + 2. * phi[6] * (r1 - c7);
    }
};

class Model8 : public simplebounce::GenericModel<8> {
  public:
    static constexpr double c0 = 0.2434;
    static constexpr double c1 = 0.5233;
    static constexpr double c2 = 0.34234;
    static constexpr double c3 = 0.4747;
    static constexpr double c4 = 0.234808;
    static constexpr double c5 = 0.57023;
    static constexpr double c6 = 0.138912;
    static constexpr double c7 = 0.51723;
    static constexpr double c8 = 0.658889;

    // potential of scalar field(s)
    double vpot(const double *phi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.) +
                     c4 * (phi[4] - 1.) * (phi[4] - 1.) +
                     c5 * (phi[5] - 1.) * (phi[5] - 1.) +
                     c6 * (phi[6] - 1.) * (phi[6] - 1.) +
                     c7 * (phi[7] - 1.) * (phi[7] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3] + phi[4] * phi[4] + phi[5] * phi[5] +
                     phi[6] * phi[6] + phi[7] * phi[7]);
        return (r1 - c8) * r2;
    }

    // derivative of potential of scalar field(s)
    void calcDvdphi(const double *phi, double *dvdphi) const {
        double r1 = (c0 * (phi[0] - 1.) * (phi[0] - 1.) +
                     c1 * (phi[1] - 1.) * (phi[1] - 1.) +
                     c2 * (phi[2] - 1.) * (phi[2] - 1.) +
                     c3 * (phi[3] - 1.) * (phi[3] - 1.) +
                     c4 * (phi[4] - 1.) * (phi[4] - 1.) +
                     c5 * (phi[5] - 1.) * (phi[5] - 1.) +
                     c6 * (phi[6] - 1.) * (phi[6] - 1.) +
                     c7 * (phi[7] - 1.) * (phi[7] - 1.));
        double r2 = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] +
                     phi[3] * phi[3] + phi[4] * phi[4] + phi[5] * phi[5] +
                     phi[6] * phi[6] + phi[7] * phi[7]);
        dvdphi[0] = 2. * c0 * (phi[0] - 1.) * r2 + 2. * phi[0] * (r1 - c8);
        dvdphi[1] = 2. * c1 * (phi[1] - 1.) * r2 + 2. * phi[1] * (r1 - c8);
        dvdphi[2] = 2. * c2 * (phi[2] - 1.) * r2 + 2. * phi[2] * (r1 - c8);
        dvdphi[3] = 2. * c3 * (phi[3] - 1.) * r2 + 2. * phi[3] * (r1 - c8);
        dvdphi[4] = 2. * c4 * (phi[4] - 1.) * r2 + 2. * phi[4] * (r1 - c8);
        dvdphi[5] = 2. * c5 * (phi[5] - 1.) * r2 + 2. * phi[5] * (r1 - c8);
        dvdphi[6] = 2. * c6 * (phi[6] - 1.) * r2 + 2. * phi[6] * (r1 - c8);
        dvdphi[7] = 2. * c7 * (phi[7] - 1.) * r2 + 2. * phi[7] * (r1 - c8);
    }
};
