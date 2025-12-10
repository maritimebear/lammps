#ifndef CRS_H
#define CRS_H

#include <vector>
#include <algorithm>
#include <cassert>

struct crs_matrix {
    std::vector<double> val;
    std::vector<size_t> col_ind;
    std::vector<size_t> row_ptr;
    std::vector<size_t> diagonal_ind; // non-zero indices of diagonal entries
    bool found_diagonal_indices = false;

    size_t nrows() const {
        return row_ptr.size() - 1; // Last element of row_ptr is nnz+1
    }

    size_t ncols() const {
        return *std::max_element(col_ind.begin(), col_ind.end()) + 1; // std::max_element() returns iterator, +1 to convert from 0-index
    }

    size_t nnz() const {
        return val.size();
    }

    void print_to_file(const std::string& filename, bool print_symmetric_entry = false) {
        // Print sparse matrix to text file
        FILE* file_handle = fopen(filename.c_str(), "w");

        // Header
        fprintf(file_handle, "%6s %6s %24s\n", "row", "col", "val");

        for (size_t row = 0; row < this->nrows(); ++row) {
            for (size_t idx_nz = row_ptr[row]; idx_nz < row_ptr[row + 1]; ++idx_nz) {
                fprintf(file_handle, "%6ld %6ld %24.15f\n", row, col_ind[idx_nz], val[idx_nz]);
                if (print_symmetric_entry) {
                    if (row != col_ind[idx_nz]) { // Avoid diagonal entries
                        fprintf(file_handle, "%6ld %6ld %24.15f\n", col_ind[idx_nz], row, val[idx_nz]);
                    }
                }
            }
        }

        fclose(file_handle);
        return;
    }

    void find_diagonal_indices() {
        if (found_diagonal_indices) {
            return;
        }
        // else find non-zero indices of diagonal entries
        diagonal_ind.reserve(nrows());
        for (size_t row = 0; row < nrows(); ++row) {
            for (size_t idx_nz = row_ptr[row]; idx_nz < row_ptr[row + 1]; ++idx_nz) {
                if (col_ind[idx_nz] == row) {
                    diagonal_ind.push_back(idx_nz); // Save non-zero index of diagonal entries
                    break;
                }
            }
        }
        assert(diagonal_ind.size() == nrows());
        found_diagonal_indices = true;
        return;
    }

};

template <typename T>
std::vector<double> crs_mvm(const crs_matrix& A, const std::vector<T>& x) {
  // CRS matrix-vector product

  std::vector<double> Ax(A.nrows(), 0.0);

  // TODO Remove after testing diagonal access
  assert(A.found_diagonal_indices);
  for (size_t row = 0; row < A.nrows(); ++row) {
      // Multiply with diagonal first
      size_t idx_nz = A.diagonal_ind[row];
      size_t col = A.col_ind[idx_nz];
      Ax[row] = (A.val[idx_nz] * x[col]);
  }
  for (size_t row = 0; row < A.nrows(); ++row) {
      for (size_t idx_nz = A.row_ptr[row]; idx_nz < A.row_ptr[row+1]; ++idx_nz) {
          size_t col = A.col_ind[idx_nz];
          if (row == col) {
              continue;
          }
          Ax[row] += (A.val[idx_nz] * x[col]);
          Ax[col] += (A.val[idx_nz] * x[row]); // Symmetric entry
      }
  }

  // for (size_t row = 0; row < A.nrows(); ++row) {
  //     for (size_t idx_nz = A.row_ptr[row]; idx_nz < A.row_ptr[row+1]; ++idx_nz) {
  //         size_t col = A.col_ind[idx_nz];
  //         Ax[row] += (A.val[idx_nz] * x[col]);
  //         if (row != col) { // Avoid multiplying diagonal twice
  //             Ax[col] += (A.val[idx_nz] * x[row]); // Symmetric entry
  //         }
  //     }
  // }

  return Ax;
}

#endif
