#include <iostream>
#include "clebsh.h"

using namespace std;

// Explicitly act with spherical part of laplace operator on spherical
// harmonic. Eigen value l*(l+1) should be returned.
complex<double> lapSphHarm(double dth, int l, int m, double theta, double phi) {
    return -(1/sin(theta)/dth*(sin(theta+dth)*(sphHarm(l, m, theta+2.*dth, phi) - sphHarm(l, m, theta+dth, phi))/dth - sin(theta)*(sphHarm(l, m, theta+dth, phi) - sphHarm(l, m, theta, phi))/dth) + 1/sin(theta)/sin(theta)*(sphHarm(l, m, theta, phi+2.*dth)+sphHarm(l, m, theta, phi)-2.*sphHarm(l, m, theta, phi+dth))/dth/dth)/sphHarm(l, m, theta, phi);
}

int main() {
    int lmax = 8;
    double dth = 1e-6;

    for (int l=0; l<lmax; l++) {
        for (int m=l; m>=-l; m--) {
            cout << l << " - " << m << " - " << lapSphHarm(dth, l, m, 0.2, 1.2) << ", " << l*(l+1) <<  endl;
        }
    }
}
