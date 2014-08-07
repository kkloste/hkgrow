/**
 * @file taydeg.hpp
 * A function to set the number of terms needed
 * in a Taylor polynomial approximation for e^1
 * to yield an error < tol.
 */

int taylordegree(const double tol) {
    double error = M_E-1;
    int N = 1;
    int k = 0;
    while (error > tol){
        k++;
        N = N*k;
        error = error - (double)1/N;
    }
    return k;
}

int alphataylordegree(const double alpha, const double tol) {
    double error = exp(alpha)-1;
    double last = 1.;
    double k = 0.;
    while(error > tol){
        k = k + 1.;
        last = (last*alpha)/k;
        error = error - last;
    }
    return (int)k;
}