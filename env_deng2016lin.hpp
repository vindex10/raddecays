class EnvLin {
public:
    // Free params
    double sigma, b, alphaS, muR, mC, rC;
    
    // Helping functions
    double smearedDelta (double r);

    // Potentials
    double Vv (double r); // Coulomb
    double Vlin(double r); // Cornell linear confinement
    double dVlin(double r); // Derivative of Vlin
    double Vss(double r, double xS, double xS1, double xS2); // Spin-spin doubleeraction
    double Vsl(double r, double xJ, double xL, double xS); // Spin-orbit doubleeraction
};
