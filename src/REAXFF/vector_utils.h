#ifndef VECUTILS_H
#define VECUTILS_H

#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>
// #include "crs_matrix.h"

template <typename T>
std::vector<size_t> sort_permutation(const std::vector<T>& vec_to_sort) {
  // Returns indices that would yield the sorted vector
  // https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
  std::vector<size_t> indices(vec_to_sort.size());
  std::iota(indices.begin(), indices.end(), 0); // Initialise vector of indices
  std::sort(indices.begin(), indices.end(), [&] (size_t i, size_t j) { return vec_to_sort[i] < vec_to_sort[j]; });
  return indices;
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
    // Add vectors elementwise
    // Not checking vector sizes
    std::vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

// template <typename T>
// std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
//     // Add vectors elementwise
//     // Not checking vector sizes
//     std::vector<T> result(v1.size());
//     std::transform(v1.cbegin(), v1.cend(), v2.cbegin(), result.begin(), std::plus<>{}); // NOT GUARANTEED TO BE IN ORDER
//     return result;
// }

template <typename T>
std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2) {
    // Subtract vectors elementwise
    // Not checking vector sizes
    std::vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

// template <typename T>
// std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2) {
//     // Subtract vectors elementwise
//     // Not checking vector sizes
//     std::vector<T> result(v1.size());
//     std::transform(v1.cbegin(), v1.cend(), v2.cbegin(), result.begin(), std::minus<>{}); // NOT GUARANTEED TO BE IN ORDER
//     return result;
// }

template <typename T>
std::vector<T> operator*(const T& scalar, const std::vector<T>& vector) {
    // Scale vector elementwise
    std::vector<T> result(vector.size());
    for (size_t i = 0; i < vector.size(); ++i) {
        result[i] = scalar * vector[i];
    }
    return result;
}

template <typename T>
double inner_product(const std::vector<T>& v1, const std::vector<T>& v2) {
    // Not checking vector sizes
    return std::inner_product(v1.cbegin(), v1.cend(), v2.cbegin(), 0.0);
    // double result = 0.0;
    // for (size_t i = 0; i < v1.size(); ++i) {
    //     result += v1[i] * v2[i];
    // }
    // return result;
}

// template <typename T>
// T operator*(const std::vector<T>& v1, const std::vector<T>& v2) {
//     // '*' operator overload for vector-vector inner product
//     return std::inner_product(v1.cbegin(), v1.cend(), v2.begin(), 0);
// }

template <typename T>
double norm(const std::vector<T>& v1) {
    return std::sqrt(inner_product(v1, v1));
}


// double operator*(const std::vector<double>& v1, const std::vector<double>& v2) {
//     // '*' operator overload for vector-vector inner product
//     return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
// }

// template <typename T>
// std::vector<double> operator*(const crs_matrix& A, const std::vector<T>& x) {
//     // '*' operator overload for csr matrix-vector product
//     return crs_mvm(A, x);
// }


#endif
