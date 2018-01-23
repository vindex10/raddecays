#include <iostream>
#include <complex>
#include <cmath>
#include "cubature.h"
#include "clebsh.h"

using namespace std;

int sphharm_f(unsigned ndim
             ,const double *x
             ,void *fdata
             ,unsigned fdim
             ,double *fval) {
    int *qnums = (int*)fdata;
    fval[0] = norm(sphHarm(qnums[0], qnums[1], x[0], x[1]))*sin(x[0]);
    return 0;
}

int main() {
    int qnums[2];
    const int maxl = 8;

    double xmin[2] = {0, 0};
    double xmax[2] = {M_PI, 2*M_PI};
    double val, err;
    for (int l=0; l<=maxl; l++) {
        for (int m=l; m>=-l; m--) {
            qnums[0] = l;
            qnums[1] = m;
            hcubature(1, sphharm_f, &qnums,  2, xmin, xmax, 0, 1e-5, 0, ERROR_INDIVIDUAL, &val, &err);
            cout << l << " - " << m << ": " << val << " ± " << err << endl;
        }
    }

    return 0;
}
