/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Implements charge-induced dipole potential for neutral-ion interactions in
// one-component plasma: Acciarri et al., "Strong Coulomb coupling influences ion
// and neutral temperatures in atmospheric pressure plasmas",
// Plasma Sources Sci. Technol. 31 (2022) 125005.
// Implementation author: John Mampilli, Imperial College London, j.mampilli@imperial.ac.uk

#include "pair_dipole_induced_acciarri.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairDipoleInducedAcciarri::PairDipoleInducedAcciarri(LAMMPS* lmp) : Pair(lmp),
    k(nullptr), r_phi_sq(nullptr), cut(nullptr), cut_global(0.0) {

    single_enable = 1;
    writedata = 1;
    respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairDipoleInducedAcciarri::~PairDipoleInducedAcciarri() {
    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(k);
        memory->destroy(r_phi_sq);
        memory->destroy(cut);
    }
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::compute(int eflag, int vflag) {
    double evdwl = 0.0; // TODO evdwl or ecoul?
    ev_init(eflag, vflag);

    double **x = atom->x;
    double *q = atom->q;
    double **f = atom->f;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    double *special_lj = force->special_lj;
    int newton_pair = force->newton_pair;

    int inum = list->inum;
    int* ilist = list->ilist;
    int* numneigh = list->numneigh;
    int** firstneigh = list->firstneigh;

    // Loop over owned atoms
    for (int ii = 0; ii < inum; ++ii) {
        int i = ilist[ii]; // Local index of current owned atom
        double xtmp = x[i][0];
        double ytmp = x[i][1];
        double ztmp = x[i][2];
        double q2tmp = q[i] * q[i];
        int itype = type[i];
        int* jlist = firstneigh[i];
        int jnum = numneigh[i];

        // Loop over neighbour atoms
        for (int jj = 0; jj < jnum; ++jj) {
            int j = jlist[jj]; // Local index of current neighbour atom
            double factor = special_lj[sbmask(j)]; // TODO special_lj or special_coul? Also change in PairDipoleInducedAcciarri::single()
            j &= NEIGHMASK;

            double delx = xtmp - x[j][0];
            double dely = ytmp - x[j][1];
            double delz = ztmp - x[j][2];
            double rsq = delx*delx + dely*dely + delz*delz;
            int jtype = type[j];

            if (rsq < cutsq[itype][jtype]) {
                double r4 = rsq * rsq; // r^4
                double r6 = r4 * rsq; // r^6
                double rratio8 = pow((r_phi_sq[itype][jtype] / rsq), 4.0); // (r_phi/r)^8
                double fpair = factor * q2tmp * (4.0 * k[itype][jtype] / r6) * (3.0 * rratio8 - 1); // Force per unit length

                f[i][0] += delx * fpair;
                f[i][1] += dely * fpair;
                f[i][2] += delz * fpair;
                if (newton_pair || j < nlocal) {
                    f[j][0] -= delx * fpair;
                    f[j][1] -= dely * fpair;
                    f[j][2] -= delz * fpair;
                }

                if (eflag) {
                    evdwl = factor * q2tmp * (k[itype][jtype] / r4) * (rratio8 - 1);
                }
                if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
            }
        }
    }
    if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::allocate() {
    // allocate all arrays
    allocated = 1;
    int np1 = atom->ntypes + 1;

    memory->create(setflag, np1, np1, "pair:setflag");
    for (int i = 1; i < np1; i++) {
        for (int j = i; j < np1; j++) {
            setflag[i][j] = 0;
        }
    }

    memory->create(cutsq, np1, np1, "pair:cutsq");
    memory->create(cut, np1, np1, "pair:cut");
    memory->create(k, np1, np1, "pair:k");
    memory->create(r_phi_sq, np1, np1, "pair:r_phi_sq");
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::settings(int narg, char** arg) {
   // global settings
   if (narg != 1) error->all(FLERR, "Pair style dipole/induced must have exactly one argument");
   cut_global = utils::numeric(FLERR, arg[0], false, lmp);

   // reset per-type pair cutoffs that have been explicitly set previously

   if (allocated) {
       for (int i = 1; i <= atom->ntypes; i++) {
           for (int j = i; j <= atom->ntypes; j++) {
               if (setflag[i][j]) cut[i][j] = cut_global;
           }
       }
   }
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::coeff(int narg, char** arg) {
    // set pair coeffs
    if (narg < 7 || narg > 8) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
    if (!allocated) allocate();

    int ilo, ihi, jlo, jhi;
    utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

    double factor = utils::numeric(FLERR, arg[2], false, lmp); // 1 / (8pi * vacuum pemittivity)
    double rel_pol = utils::numeric(FLERR, arg[3], false, lmp); // relative polarisability
    double bohr_radius = utils::numeric(FLERR, arg[4], false, lmp);
    // Coefficients to calculate r_phi = c * a_in
    double c = utils::numeric(FLERR, arg[5], false, lmp);
    double a_in = utils::numeric(FLERR, arg[6], false, lmp);
    double cut_tmp = cut_global;
    if (narg == 10) cut_tmp = utils::numeric(FLERR, arg[7], false, lmp);

    // Calculate k and r_phi_sq
    double k_tmp = factor * rel_pol * pow(bohr_radius, 3.0);
    double r_phi_sq_tmp = pow((c * a_in), 2.0);

    int count = 0;
    for (int i = ilo; i <= ihi; ++i) {
        for (int j = MAX(jlo, i); j <= jhi; ++j) {
            cut[i][j] = cut_tmp;
            k[i][j] = k_tmp;
            r_phi_sq[i][j] = r_phi_sq_tmp;
            setflag[i][j] = 1;
            ++count;
        }
    }

    if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ---------------------------------------------------------------------- */

double PairDipoleInducedAcciarri::init_one(int i, int j) {
    // init for one type pair i,j and corresponding j,i
    if (setflag[i][j] == 0) {
        error->all(FLERR, Error::NOLASTLINE,
                   "All pair coeffs are not set. Status:\n" + Info::get_pair_coeff_status(lmp));
    }

    if (offset_flag) {
        error->all(FLERR, Error::NOLASTLINE,
                "pair_modify shift yes not implemented for pair_style dipole/induced");
    }

    k[j][i] = k[i][j];
    r_phi_sq[j][i] = r_phi_sq[i][j];

    return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::write_restart(FILE* fp) {
    // proc 0 writes to restart file
    write_restart_settings(fp);

    for (int i = i; i <= atom->ntypes; ++i) {
        for (int j = i; j <= atom->ntypes; ++j) {
            fwrite(&setflag[i][j], sizeof(int), 1, fp);
            if (setflag[i][j]) {
                fwrite(&k[i][j], sizeof(double), 1, fp);
                fwrite(&r_phi_sq[i][j], sizeof(double), 1, fp);
                fwrite(&cut[i][j], sizeof(double), 1, fp);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::read_restart(FILE* fp) {
    // proc 0 reads from restart file, bcasts
    read_restart_settings(fp);
    allocate();
    int me = comm->me;
    for (int i = 1; i <= atom->ntypes; ++i) {
        for (int j = i; j <= atom->ntypes; ++j) {
            if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
            MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
            if (setflag[i][j]) {
                if (me == 0) {
                    utils::sfread(FLERR, &k[i][j], sizeof(double), 1, fp, nullptr, error);
                    utils::sfread(FLERR, &r_phi_sq[i][j], sizeof(double), 1, fp, nullptr, error);
                    utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
                }
                MPI_Bcast(&k[i][j], 1, MPI_DOUBLE, 0, world);
                MPI_Bcast(&r_phi_sq[i][j], 1, MPI_DOUBLE, 0, world);
                MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::write_restart_settings(FILE* fp) {
    // proc 0 writes to restart file
    fwrite(&cut_global, sizeof(double), 1, fp);
    fwrite(&offset_flag, sizeof(int), 1, fp);
    fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::read_restart_settings(FILE* fp) {
    // proc 0 reads from restart file, bcasts
    if (comm->me == 0) {
        utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
        utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
        utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    }
    MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::write_data(FILE* fp) {
    // proc 0 writes to data file
    for (int i = 1; i <= atom->ntypes; ++i) {
        fprintf(fp, "%d %g %g\n", i, k[i][i], r_phi_sq[i][i]);
    }
}

/* ---------------------------------------------------------------------- */

void PairDipoleInducedAcciarri::write_data_all(FILE* fp) {
    // proc 0 writes all pairs to data file
    for (int i = 1; i <= atom->ntypes; ++i) {
        for (int j = i; j <= atom->ntypes; ++j) {
            fprintf(fp, "%d %d %g %g %g\n", i, j, k[i][j], r_phi_sq[i][j], cut[i][j]);
        }
    }
}

/* ---------------------------------------------------------------------- */

double PairDipoleInducedAcciarri::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                                 double /*factor_coul*/, double factor_lj, double &fforce) {
    // Compute force and energy for a single pair of atoms
    // Force per unit distance written to double& fforce, energy returned

    double r4 = rsq * rsq; // r^4
    double r6 = r4 * rsq; // r^6
    double rratio8 = pow((r_phi_sq[itype][jtype] / rsq), 4.0); // (r_phi/r)^8
    double factor = factor_lj; // TODO: factor_lj or factor_coul?

    fforce = factor * (4.0 * k[itype][jtype] / r6) * (3.0 * rratio8 - 1); // Force per unit length
    return factor * (k[itype][jtype] / r4) * (rratio8 - 1); // Energy
}

/* ---------------------------------------------------------------------- */

void* PairDipoleInducedAcciarri::extract(const char* str, int& dim) {
    dim = 2; // what is this
    if (strcmp(str, "k") == 0) return (void *) k;
    if (strcmp(str, "r_phi_sq") == 0) return (void *) r_phi_sq;

    return nullptr;
}
