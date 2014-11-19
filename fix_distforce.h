#ifdef FIX_CLASS

FixStyle(distforce,FixDistForce)

#else

#ifndef LMP_FIX_DISTFORCE_H
#define LMP_FIX_DISTFORCE_H

#include "fix.h"
#include "edm_bias.h"
#include "random_mars.h"
#include "compute.h"

namespace LAMMPS_NS {

class FixDistForce : public Fix {
 public:
  FixDistForce(class LAMMPS *, int, char **);
  ~FixDistForce();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  void post_integrate();
  void min_pre_exchange();
//  void end_of_step();

 private:

 int nmolecules;
 tagint idlo,idhi;
 int firstflag;
// double *counterproc,*counterall;
 double *massproc,*masstotal;
 double **comhold, **comall;
 double **forceproc, **forceall;
// void summolmass(int);

};

};

#endif
#endif

