/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sw/wfnho,PairSWWFNHO)

#else

#ifndef LMP_PAIR_SW_WFNHO_H
#define LMP_PAIR_SW_WFNHO_H

#include "pair.h"
#define PIVAL 3.1415926535898

namespace LAMMPS_NS {

class PairSWWFNHO : public Pair {
 public:
  PairSWWFNHO(class LAMMPS *);
  ~PairSWWFNHO();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  
 private:
  struct Param {
    double epsilon,sigma;
    double aij,aik,cij,lambda,gammaij,gammaik,costheta;
    double biga,bigb;
    double powerp,powerq;
    double tol;
    double cutpair,cutij,cutik,cutpairsq,cutijsq,cutiksq;
    double sigma_gammaij,sigma_gammaik,lambda_epsilon,lambda_epsilon2;
    double c1,c2,c3,c4,c5,c6;
    int ielement,jelement,kelement;
  };
   struct Softparam {
    double ma,mb,mc,md,me,bigr,bigd;
    int ielement,jelement;
  };
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int **elem2soft;   		  // mapping from element doubles to bond softening
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction
  int nsofts;			  // # of bond softening parameter sets
  int maxsofts;			  // max # of softening parameter sets
  Softparam *soft;              // parameter set for the bond softening
  int **softflag;
  double *coord;
  void allocate();
  void read_file(char *);
  void setup();
  void twobody(Param *, double, double &, int, double &);
  void threebody(Param *, double, double, double *, double *,
		 double *, double *, int, double &);
  double gsoft(double, int, int);
};

}

#endif
#endif