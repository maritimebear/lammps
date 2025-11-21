#ifndef LINEAR_SOLVERS_H
#define LINEAR_SOLVERS_H

#include <vector>
#include "crs_matrix.h"
#include "vector_utils.h"


template <typename DUMMY>
int CRS_BiCGStab(const crs_matrix& A, std::vector<double>& x, const std::vector<double>& b, double tolerance, double rhotol, int maxiters) {
    // BiCGStab using CRS matrix and std::vectors
    // Returns iteration count
    // Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods, Figure 2.10
    
    DUMMY dummy; // To keep implementation in .h file

    printf("A.nrows(): %ld\n", A.nrows());
    printf("x.size(): %ld\n", x.size());
    printf("b.size(): %ld\n", b.size());

    double bnorm2 = std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
    if (bnorm2 == 0.0) {
        // error->warning(FLERR, "BiCGStab(): ||b|| == 0.0, b == zero vector?");
        printf("BiCGStab(): ||b|| == 0.0, b == zero vector?\n");
        return 0;
    }

    std::vector<double> r = b - crs_mvm(A, x);
    printf("r.size(): %ld\n", r.size());

    // Convergence check
    double rnorm2 = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
    if (rnorm2 < bnorm2 * (tolerance * tolerance)) {
        return 0;
    }

    std::vector<double> rhat = r; // Shadow residual
    printf("rhat.size(): %ld\n", rhat.size());

    // Variables declared here since referenced in loop before assignment
    std::vector<double> p(x.size());
    std::vector<double> v(x.size());
    double rho_old, alpha, omega;

    for (int iter = 1; iter < maxiters; ++iter) {

        double rho = std::inner_product(rhat.begin(), rhat.end(), r.begin(), 0.0);
        if (fabs(rho) < rhotol) {
            // error->warning(FLERR, "BiCGStab(): |rho| = {:.2} < rhotol = {:.2}", fabs(rho), rhotol);
            printf("BiCGStab(): |rho| < rhotol\n");
            break;
        }

        if (iter == 1) {
            p = r;
        } else {
            double beta = (rho / rho_old) * (alpha / omega);
            p = r + (beta * (p - (omega * v)));
        }

        // TODO preconditioning: p_hat = M^-1 * p
        std::vector<double> p_hat = p;

        v = crs_mvm(A, p_hat);
        printf("v.size(): %ld\n", v.size());

        double rhat_v = std::inner_product(rhat.begin(), rhat.end(), v.begin(), 0.0);
        if (rhat_v < rhotol) {
            // error->warning(FLERR, "BiCGStab(): <rhat, v> = {:.2} < rhotol = {:.2}", rhat_v, rhotol);
            printf("BiCGStab(): <rhat, v> < rhotol\n");
            break;
        }

        alpha = rho / rhat_v;
        printf("alpha: %f\n", alpha);

        std::vector<double> s = r - (alpha * v);
        printf("s.size(): %ld\n", s.size());

        // Convergence check
        if (std::inner_product(s.begin(), s.end(), s.begin(), 0.0) < (tolerance * tolerance)) {
            x = x + (alpha * p_hat);
            return iter;
        }

        // TODO preconditioning: s_hat = M^1 * s
        std::vector<double> s_hat = s;
        printf("s_hat.size(): %ld\n", s_hat.size());

        std::vector<double> t = crs_mvm(A, s_hat);

        omega = std::inner_product(t.begin(), t.end(), s.begin(), 0.0) / (t.begin(), t.end(), t.begin(), 0.0);
        if (fabs(omega) < rhotol) {
            // error->warning(FLERR, "BiCGStab(): |omega| = {:.2} < rhotol = {:.2}", fabs(omega), rhotol);
            printf("BiCGStab(): |omega| < rhotol\n");
            break;
        }

        x = x + ((alpha * p_hat) + (omega * s_hat));

        r = s - (omega * t);

        rnorm2 = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
        if (rnorm2 < bnorm2 * (tolerance * tolerance)) {
            return 0;
        }

        rho_old = rho;
    }

    return -1; // Only in case of numerical breakdown inside iteration loop
}

#endif
