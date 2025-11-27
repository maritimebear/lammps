/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(acks2/reax,FixACKS2ReaxFF);
FixStyle(acks2/reaxff,FixACKS2ReaxFF);
// clang-format on
#else

#ifndef LMP_FIX_ACKS2_REAXFF_H
#define LMP_FIX_ACKS2_REAXFF_H

#include "fix_qeq_reaxff.h"
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include "crs_matrix.h"
#include "vector_utils.h"

namespace LAMMPS_NS {

class FixACKS2ReaxFF : public FixQEqReaxFF {
 public:
  FixACKS2ReaxFF(class LAMMPS *, int, char **);
  ~FixACKS2ReaxFF() override;
  void post_constructor() override;
  void init() override;
  void init_storage() override;
  void pre_force(int) override;

  double *get_s() { return s; }

  double compute_scalar() override;

 protected:
  int NN, last_rows_rank, last_rows_flag;

  double **s_hist_X, **s_hist_last;
  double *bcut_acks2, bond_softness, **bcut;    // acks2 parameters
  // bcut_acks2: atom-type specific cutoff parameter, bond_softness: global bond softness parameter, bcut: average of 2-atom type bcut_acks2 parameters

  sparse_matrix X;
  double *Xdia_inv;
  double *X_diag;

  //BiCGStab storage
  double *g, *q_hat, *r_hat, *y, *z;

  // TODO remove test/debug variables
  bool print_system;
  bool print_acks2_matrix;
  std::vector<double> vec_b_s;
  std::vector<double> vec_s;
  std::vector<double> vec_H_diag;
  std::vector<double> vec_Hdia_inv;
  std::vector<double> vec_X_diag;
  std::vector<double> vec_Xdia_inv;

  void pertype_parameters(char *) override;
  void init_bondcut();
  void allocate_storage() override;
  void deallocate_storage() override;
  void allocate_matrix() override;
  void deallocate_matrix() override;

  void init_matvec() override;
  void compute_X();
  // bool Xflag(int, int); // TODO Cleanup
  double calculate_X(double, double) const;
  void calculate_Q() override;

  int BiCGStab(double *, double *);
  void sparse_matvec_acks2(sparse_matrix *, sparse_matrix *, double *, double *);

  // TODO Cleanup
  std::string append_timestep(const std::string&);
  void print_sparse_matrix(sparse_matrix&, const std::string&); // TODO: Remove after debugging ACKS2?
  void print_matrix_diagonals();
  void print_array(double*, int, const std::string&);
  // int copy_array_to_vector(double* , std::vector<double>&);
  void copy_X_diag(std::vector<double>&);
  void copy_Xdia_inv(std::vector<double>&);
  int copy_H_diag(std::vector<double>&);
  int copy_Hdia_inv(std::vector<double>&);

  int ACKS2CG(double*, double*);
  int ACKS2BiCGStab(double*, double*, double, int);
  int BiCGStab_NoComm(double*, double*, double, int);
  int _ACKS2BiCGStab(double*, double*, double, int);
  int CRS_BiCGStab(const crs_matrix&, std::vector<double>&, const std::vector<double>&, double, double, int) const;
  int CRS_CG(const crs_matrix&, std::vector<double>&, const std::vector<double>&, double, double, int) const;
  int RestartedBiCGStab(double*, double*, double, int);
  void copy_array_to_vector(double*, std::vector<double>&) const;
  void compare_vectors(std::vector<double>&, std::vector<double>&) const;
  std::vector<double> construct_acks2_rhs(double*) const;
  std::vector<double> array_to_vector(double*) const;
  void vector_to_array(std::vector<double>&, double*, const std::unordered_map<int, int>&);
  void compare_array_vector(double*, const std::vector<double>&, double) const;
  std::unordered_map<int, int> construct_tag_map() const;
  crs_matrix assemble_acks2_matrix(const std::unordered_map<int, int>&) const;
  // crs_matrix __assemble_acks2_matrix(const std::unordered_map<int, int>&) const;

  bool array_vec_equal(double*, const std::vector<double>&);
  bool diag_vec_equal(double*, const std::vector<double>&);

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void more_forward_comm(double *);
  void more_reverse_comm(double *);
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double parallel_norm(double *, int) override;
  double parallel_dot(double *, double *, int) override;
  double parallel_vector_acc(double *, int) override;

  void vector_sum(double *, double, double *, double, double *, int) override;
  void vector_add(double *, double, double *, int) override;
  void vector_copy(double *, double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
