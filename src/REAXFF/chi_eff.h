/* Hardcoding Gaussian orbital exponents for calculating effective electronegativities in ACKS2 with electric fields and periodic boundaries
 * Following QEqR/QTPIE
 * Orbital exponents from Table 2.2, Jiahao Chen, "Theory and Applications of Fluctuating-Charge Models", 2010, arXiv:1004.0186 
*/

#ifndef CHI_EFF
#define CHI_EFF

#include <map>
#include <string>
#include <algorithm>
#include <cmath>
#include <utility>

namespace chi_eff {

constexpr double ANGSTROM_TO_BOHRRADIUS_SQ = 3.571064831;
constexpr double scale = 1.0;

const std::map<std::string, double> gauss_exp {
    // Table 2.2, Jiahao Chen, "Theory and Applications of Fluctuating-Charge Models", 2010, arXiv:1004.0186 
    {"H", 0.5434 * ANGSTROM_TO_BOHRRADIUS_SQ},
    {"O", 0.2240 * ANGSTROM_TO_BOHRRADIUS_SQ}
};


// Calculate a cutoff distance to neglect overlap integrals in calc_chi_eff() when less than pow(10, -olap_cut)

constexpr int olap_cut = 10;

const double exp_min = std::min_element(gauss_exp.cbegin(), gauss_exp.cend(),
        [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; }
        )->second; // std::min_element() -> std::map::iterator, iterator to pair of <key, value>

const double dist_cutoff_sq = 2 * olap_cut * std::log(10.0) / exp_min;


// Functions to calculate terms in overlap integral

constexpr double calculate_expfactor(double exp_i, double exp_j) {
    return exp_i * exp_j / (exp_i + exp_j); 
}

constexpr double calculate_prefactor(double exp_i, double exp_j) {
    return std::pow((4.0 * calculate_expfactor(exp_i, exp_j) / (exp_i + exp_j)), 0.75);
}


}

#endif
