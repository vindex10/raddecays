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
    complex<double> prod = conj(sphHarm(qnums[0], qnums[1], x[0], x[1]))*sphHarm(qnums[2], qnums[3], x[0], x[1])*sin(x[0]);
    fval[0] = real(prod);
    fval[1] = imag(prod);
    return 0;
}

int main() {
    int qnums[4];
    const int maxl = 8;

    double xmin[2] = {0, 0};
    double xmax[2] = {M_PI, 2*M_PI};
    double val[2], err;
    for (int l1=0; l1<=maxl; l1++) {
        for (int m1=l1; m1>=-l1; m1--) {
            for (int l2=0; l2 <=l1; l2++) {
                for (int m2=l2; m2 >= -l2; m2--) {
                    qnums[0] = l1;
                    qnums[1] = m1;
                    qnums[2] = l2;
                    qnums[3] = m2;
                    hcubature(2, sphharm_f, &qnums,  2, xmin, xmax, 0, 1e-5, 0, ERROR_INDIVIDUAL, val, &err);
                    cout << "(" << l1 << "," << m1 << "), (" << l2 << "," << m2 << ")" << ": " << val[0] << "+ i" << val[1] << " ± " << err << endl;
                }
            }
        }
    }

    return 0;
}
