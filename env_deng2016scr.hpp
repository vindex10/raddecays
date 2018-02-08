class EnvScr {
public:
    // Free params
    double sigma, b, alphaS, muR, mC, mu, rC;
    
    // Helping functions
    double smearedDelta (double r);

    // Potentials
    double Vv (double r); // Coulomb
    double Vscr(double r); // Exponentially screened confinement
    double dVscr(double r); // Derivative of Vscr
    double Vss(double r, double xS, double xS1, double xS2); // Spin-spin interaction
    double Vsl(double r, double xJ, double xL, double xS); // Spin-orbit interaction
};
