
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation. Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software. This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Lucas Hale
   Modified from the Stillinger-Weber potential by: Aidan Thompson (SNL)
   
   Note 2015/09/16
   This pair style and coefficient data set (total 3 files) was extracted from LM Hale phd thesis;
   this text-compilable version was updated on 2015/09/16 for lammps version 20150810 by Alessandro L. Sellerio;
   the files may contain errors and typos introduced by me due to pdf to text conversion procedure.
   All credits go to original authors.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_sw_wfnho.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairSWWFNHO::PairSWWFNHO(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;
  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  nsofts = 0;
  params = NULL;
  elem2param = NULL;
  soft = NULL;
  elem2soft = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairSWWFNHO::~PairSWWFNHO()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(soft);
  memory->destroy(elem2param);
  
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(softflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] coord;
  }
}

/* ---------------------------------------------------------------------- */
void PairSWWFNHO::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1,itag,jtag;
  int itype,jtype,ktype,ijparam,ijkparam;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2,bigr,bigd,gij,r;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  double **x = atom->x;
  double **f = atom->f;
  int *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // calculate coordination number for softening function

for (ii = 0; ii < inum; ii++)
{
  i = ilist[ii];
  coord[tag[i]]=0;
}

for (ii = 0; ii < inum; ii++) {
  i = ilist[ii];
  itag = tag[i];
  itype = map[type[i]];
  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];
  jlist = firstneigh[i];
  jnum = numneigh[i];
  for (jj = 0; jj < jnum; jj++) {
    j=jlist[jj];
    jtag = tag[j];
    jtype = map[type[j]];
    if (softflag[itype][jtype])
          {
            ijparam = elem2param[itype][jtype][jtype];
            bigr = soft[elem2soft[itype][jtype]].bigr;
            bigd = soft[elem2soft[itype][jtype]].bigd;
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            r = sqrt(rsq)/params[ijparam].sigma;
            if (r < (bigr - bigd)) coord[itag] += 1;
            else if (r < (bigr + bigd))
              coord[itag] += 1 - (r-bigr+bigd)/(2*bigd) + sin(PIVAL*(r-bigr+bigd)/bigd)/(2*PIVAL);
          }
  }
}

  // loop over full neighbor list of my atoms

for (ii = 0; ii < inum; ii++) {
  i = ilist[ii];
  itag = tag[i];
  itype = map[type[i]];
  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];

  // two-body interactions, skip half of them

  jlist = firstneigh[i];
  jnum = numneigh[i];

  for (jj = 0; jj < jnum; jj++) {
  j = jlist[jj];
  jtag = tag[j];

  if (itag > jtag) {
         if ((itag+jtag) % 2 == 0) continue;
  } else if (itag < jtag) {
         if ((itag+jtag) % 2 == 1) continue;
  } else {
         if (x[j][2] < ztmp) continue;
         if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
         if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
  }

  jtype = map[type[j]];
 
  delx = xtmp - x[j][0];
  dely = ytmp - x[j][1];
  delz = ztmp - x[j][2];
  rsq = delx*delx + dely*dely + delz*delz;
  
  ijparam = elem2param[itype][jtype][jtype];
  if (rsq > params[ijparam].cutpairsq) continue;
  
  twobody(&params[ijparam],rsq,fpair,eflag,evdwl);
  
  if (softflag[itype][jtype]) gij = gsoft(coord[itag],itype,jtype);
  else if (softflag[jtype][itype]) gij = gsoft(coord[jtag],jtype,itype);
  else gij = 1;
  
  evdwl = gij*evdwl;
  fpair = gij*fpair;
  
  f[i][0] += delx*fpair;
  f[i][1] += dely*fpair;
  f[i][2] += delz*fpair;
  f[j][0] -= delx*fpair;
  f[j][1] -= dely*fpair;
  f[j][2] -= delz*fpair;
  if (evflag) ev_tally(i,j,nlocal,newton_pair,
                               evdwl,0.0,fpair,delx,dely,delz);
  }

  jnumm1 = jnum - 1;

  for (jj = 0; jj < jnumm1; jj++) {
  j = jlist[jj];
  jtype = map[type[j]];
  delr1[0] = x[j][0] - xtmp;
  delr1[1] = x[j][1] - ytmp;
  delr1[2] = x[j][2] - ztmp;
  rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
  
  for (kk = jj+1; kk < jnum; kk++) {
    k = jlist[kk];
    ktype = map[type[k]];
    ijkparam = elem2param[itype][jtype][ktype];
    
    if (rsq1 > params[ijkparam].cutijsq) continue;
    
    delr2[0] = x[k][0] - xtmp;
    delr2[1] = x[k][1] - ytmp;
    delr2[2] = x[k][2] - ztmp;
    rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
    
    if (rsq2 > params[ijkparam].cutiksq) continue;
    
    threebody(&params[ijkparam],rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);
    
    f[i][0] -= fj[0] + fk[0];
    f[i][1] -= fj[1] + fk[1];
    f[i][2] -= fj[2] + fk[2];
    f[j][0] += fj[0];
    f[j][1] += fj[1];
    f[j][2] += fj[2];
    f[k][0] += fk[0];
    f[k][1] += fk[1];
    f[k][2] += fk[2];
    if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
    }
   }
  }
  
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

double PairSWWFNHO::gsoft(double cn, int i, int j)
{
           double ma,mb,mc,md,me,first,second;
           int ij = elem2soft[i][j];
           ma = soft[ij].ma;
           mb = soft[ij].mb;
           mc = soft[ij].mc;
           md = soft[ij].md;
           me = soft[ij].me;
           first = ma/(exp((mb-cn)/mc)+1);
           second = exp(md*(cn-me)*(cn-me));
           return first*second;
}

/* ---------------------------------------------------------------------- */
void PairSWWFNHO::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  
  int num = ceil(atom->natoms);
  coord = new double[num+1];
  
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(softflag,n+1,n+1,"pair:softflag");
  
  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSWWFNHO::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSWWFNHO::coeff(int narg, char **arg)
{
  int i,j,n;
  
  if (!allocated) allocate();
  
  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters
  
  read_file(arg[2]);
  setup();
  
  // clear setflag since coeff() called once with I,J = * *
  
  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
  
  // set setflag i,j for type pairs where both are mapped to elements
  
  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
            setflag[i][j] = 1;
            count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSWWFNHO::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style Stillinger-Weber requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style Stillinger-Weber requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSWWFNHO::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairSWWFNHO::read_file(char *file)
{
  int params_per_line = 17;
  int soft_per_line = 10;
  int per_line;
  int param_type = 0;
  char **words = new char*[params_per_line+1];
  char look1,look2,look3;
  
  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;
  
  if (soft != NULL) free(soft);
  soft = NULL;
  nsofts = maxsofts = 0;
  for (int y=0;y<atom->ntypes;y++)
    for (int z=0;z<atom->ntypes;z++)
      softflag[y][z]=0;
  
  // open file on proc 0
  
  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open Stillinger-Weber potential file %s",file);
      error->one(FLERR,str);
    }
  }
  
  // read each set of params from potential file
  // one set of params can span multiple lines
  // store params if all 3 element tags are in element list
  
  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;
  
  while (1) {
    param_type=0;
    per_line = params_per_line;
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
              eof = 1;
              fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank
    if (ptr = strchr(line,'#')) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // determine which parameter type line consists of
    if ((ptr=strchr(line,'S'))) {
      look1 = *(ptr+1);
      look2 = *(ptr+2);
      look3 = *(ptr+3);
      if (look1 == 'o' && look2 == 'f' && look3 == 't') {
        per_line = soft_per_line;
        param_type = 1;
      }
    }
    
    while (nwords < per_line) {
    n = strlen(line);
    if (comm->me == 0) {
      ptr = fgets(&line[n],MAXLINE-n,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
     }
     MPI_Bcast(&eof,1,MPI_INT,0,world);
     if (eof) break;
     MPI_Bcast(&n,1,MPI_INT,0,world);
     MPI_Bcast(line,n,MPI_CHAR,0,world);
     if (ptr = strchr(line,'#')) *ptr = '\0';
     nwords = atom->count_words(line);
   }

   if (nwords != per_line) 
     error->all(FLERR,"Incorrect format in Stillinger-Weber potential file");
  
   // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while (words[nwords++] = strtok(NULL," \t\n\r\f")) continue;
    
    //2 and 3 body parameter settings
    if (param_type == 0) {
    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next entry in file

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;
    
    // load up parameter settings and error check their values

    if (nparams == maxparam) {
    maxparam += DELTA;
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                               "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].epsilon = atof(words[3]);
    params[nparams].sigma = atof(words[4]);
    params[nparams].aij = atof(words[5]);
    params[nparams].aik = atof(words[6]);
    params[nparams].lambda = atof(words[7]);
    params[nparams].gammaij = atof(words[8]);
    params[nparams].gammaik = atof(words[9]);
    params[nparams].costheta = atof(words[10]);
    params[nparams].biga = atof(words[11]);
    params[nparams].bigb = atof(words[12]);
    params[nparams].powerp = atof(words[13]);
    params[nparams].powerq = atof(words[14]);
    params[nparams].cij = atof(words[15]);
    params[nparams].tol = atof(words[16]);
/*
 if (params[nparams].epsilon < 0.0 || params[nparams].sigma < 0.0 ||
          params[nparams].aij < 0.0 || params[nparams].aik < 0.0 ||
          params[nparams].lambda < 0.0 || params[nparams].gammaij < 0.0 ||
          params[nparams].gammaik < 0.0 || ||
          params[nparams].bigb < 0.0 || params[nparams].powerp < 0.0 ||
          params[nparams].powerq < 0.0 || params[nparams].cij < 0.0 ||
          params[nparams].tol < 0.0)
   error->all("Illegal Stillinger-Weber parameter");
*/

    nparams++;
  }
  
  else if (param_type == 1) {
    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[1],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[2],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    
    // load up parameter settings and error check their values
    softflag[ielement][jelement] = 1;
    if (nsofts == maxsofts) {
      maxsofts += DELTA;
      soft = (Softparam *) memory->srealloc(soft,maxsofts*sizeof(Softparam),
                                                        "pair:soft");
    }
    
    soft[nsofts].ielement = ielement;
    soft[nsofts].jelement = jelement;
    soft[nsofts].ma = atof(words[3]);
    soft[nsofts].mb = atof(words[4]);
    soft[nsofts].mc = atof(words[5]);
    soft[nsofts].md = atof(words[6]);
    soft[nsofts].me = atof(words[7]);
    soft[nsofts].bigr = atof(words[8]);
    soft[nsofts].bigd = atof(words[9]);
    nsofts++;
  }
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */
void PairSWWFNHO::setup()
{
  int i,j,k,m,n,o,p;
  double rtmp1,rtmp2,rtmp3;

  // set elem2param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  if (elem2param) memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");
  
  if (elem2soft) memory->destroy(elem2soft);
  memory->create(elem2soft,nelements,nelements,"pair:elem2soft");
  
  for (i = 0; i < nelements; i++) {
    for (j = 0; j < nelements; j++) {
      o = -1;
      for (p = 0; p < nsofts; p++) {
            if (i == soft[p].ielement && j == soft[p].jelement) {
              if (o >= 0) error->all(FLERR,"Potential file has duplicate entry");
              o = p;
            }
      }
      elem2soft[i][j] = o;
      for (k = 0; k < nelements; k++) {
            n = -1;
            for (m = 0; m < nparams; m++) {
              if (i == params[m].ielement && j == params[m].jelement &&
                  k == params[m].kelement) {
                if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
                n = m;
           }
         }
         if (n < 0) error->all(FLERR,"Potential file is missing an entry");
         elem2param[i][j][k] = n;
      }
    }
  }

  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cuts must remain unchanged as it is a potential parameter
  // (cut = a*sigma)

  for (m = 0; m < nparams; m++) {
    params[m].cutpair = params[m].sigma*params[m].cij;
    params[m].cutij = params[m].sigma*params[m].aij;
    params[m].cutik = params[m].sigma*params[m].aik;
 
    rtmp1 = params[m].cutpair;
    rtmp2 = params[m].cutij;
    rtmp3 = params[m].cutik;
    if (params[m].tol > 0.0) error->all(FLERR,"Potential not currently set to accept tol values");
    
    params[m].cutpairsq = rtmp1 * rtmp1;
    params[m].cutijsq = rtmp2 * rtmp2;
    params[m].cutiksq = rtmp3 * rtmp3;
    params[m].sigma_gammaij = params[m].sigma*params[m].gammaij;
    params[m].sigma_gammaik = params[m].sigma*params[m].gammaik;
    params[m].lambda_epsilon = params[m].lambda*params[m].epsilon;
    params[m].lambda_epsilon2 = 2.0*params[m].lambda*params[m].epsilon;
    params[m].c1 = params[m].biga*params[m].epsilon *
      params[m].powerp*params[m].bigb *
      pow(params[m].sigma,params[m].powerp);
    params[m].c2 = params[m].biga*params[m].epsilon*params[m].powerq *
     pow(params[m].sigma,params[m].powerq);
    params[m].c3 = params[m].biga*params[m].epsilon*params[m].bigb *
      pow(params[m].sigma,params[m].powerp+1.0);
    params[m].c4 = params[m].biga*params[m].epsilon *
      pow(params[m].sigma,params[m].powerq+1.0);
    params[m].c5 = params[m].biga*params[m].epsilon*params[m].bigb *
      pow(params[m].sigma,params[m].powerp);
    params[m].c6 = params[m].biga*params[m].epsilon *
      pow(params[m].sigma,params[m].powerq);
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp1 = sqrt(params[m].cutpairsq);
    if (rtmp1 > cutmax) cutmax = rtmp1;
    rtmp2 = sqrt(params[m].cutijsq);
    if (rtmp2 > cutmax) cutmax = rtmp2;
    rtmp3 = sqrt(params[m].cutiksq);
    if (rtmp3 > cutmax) cutmax = rtmp3;
  }
}

/* ---------------------------------------------------------------------- */
void PairSWWFNHO::twobody(Param *param, double rsq, double &fforce,
                         int eflag, double &eng)
{
  double r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;
  
  r = sqrt(rsq);
  rinvsq = 1.0/rsq;
  rp = pow(r,-param->powerp);
  rq = pow(r,-param->powerq);
  rainv = 1.0 / (r - param->cutpair);
  rainvsq = rainv*rainv*r;
  expsrainv = exp(param->sigma * rainv);
  fforce = (param->c1*rp - param->c2*rq +
              (param->c3*rp -param->c4*rq) * rainvsq) * expsrainv * rinvsq;
  if (eflag) eng = (param->c5*rp - param->c6*rq) * expsrainv;
}

/* ---------------------------------------------------------------------- */

void PairSWWFNHO::threebody(Param *paramijk,
                           double rsq1, double rsq2,
                           double *delr1, double *delr2,
                           double *fj, double *fk, int eflag, double &eng)
{
  double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  double rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2;
  double facang,facang12,csfacang,csfac1,csfac2;
  
  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - paramijk->cutij);
  gsrainv1 = paramijk->sigma_gammaij * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = exp(gsrainv1);
  
  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  rainv2 = 1.0/(r2 - paramijk->cutik);
  gsrainv2 = paramijk->sigma_gammaik * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2;
  expgsrainv2 = exp(gsrainv2);
  
  rinv12 = 1.0/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  delcs = cs - paramijk->costheta;
  delcssq = delcs*delcs;
  
  facexp = expgsrainv1*expgsrainv2;
  
  // facrad = sqrt(paramijk->lambda_epsilon*paramijk->lambda_epsilon) *
  //       facexp*delcssq;
  
  facrad = paramijk->lambda_epsilon * facexp*delcssq;
  frad1 = facrad*gsrainvsq1;
  frad2 = facrad*gsrainvsq2;
  facang = paramijk->lambda_epsilon2 * facexp*delcs;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;
  
  fj[0] = delr1[0]*(frad1+csfac1)-delr2[0]*facang12;
  fj[1] = delr1[1]*(frad1+csfac1)-delr2[1]*facang12;
  fj[2] = delr1[2]*(frad1+csfac1)-delr2[2]*facang12;
  
  csfac2 = rinvsq2*csfacang;
  
  fk[0] = delr2[0]*(frad2+csfac2)-delr1[0]*facang12;
  fk[1] = delr2[1]*(frad2+csfac2)-delr1[1]*facang12;
  fk[2] = delr2[2]*(frad2+csfac2)-delr1[2]*facang12;
  
  if (eflag) eng = facrad;
}
