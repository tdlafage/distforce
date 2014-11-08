// fix fix_id group_id(this should be the virtual particles) distforce

#include "fix_distforce.h"
#include "error.h"
#include "atom.h"
#include "memory.h"
#include "group.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixDistForce::FixDistForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
   if (narg != 2) error->all(FLERR,"Illegal fix distforce command");

   if (atom->molecular == 0)
    error->all(FLERR,"Fix distforce requires full or molecular atom style");


  // memory->create(masstotal,nmolecules,"distforce/molecule:masstotal");

}

/* ---------------------------------------------------------------------- */

FixDistForce::~FixDistForce()
{
  memory->destroy(masstotal);
}

/* ---------------------------------------------------------------------- */

using namespace FixConst;
int FixDistForce::setmask()
{
  // set with a bitmask how and when apply the force from EDM 
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDistForce::init()
{

}

/* ---------------------------------------------------------------------- */

void FixDistForce::setup(int vflag)
{
}

/* ---------------------------------------------------------------------- */

void FixDistForce::min_setup(int vflag)
{
  setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixDistForce::pre_force(int vflag)
{
}
/* ---------------------------------------------------------------------- */

void FixDistForce::post_force(int vflag)
{

}

/* ---------------------------------------------------------------------- */

void FixDistForce::min_pre_force(int vflag)
{
 pre_force(vflag);

  tagint imol;
  double massone;
  int j;
  double masstotal;
  double unwrap[3];

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;


  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      x[i][0] = x[i][1] = x[i][2] = 0.0;
      imol = molecule[i];
//    if (molmap) imol = molmap[imol-idlo];
//    else imol--;
//    domain->unmap(x[i],image[i],unwrap);
//    if (rmass) massone = rmass[i];
//    massone = mass[type[i]];
      masstotal = 0.0;
      for(int j = 0; j < nlocal; j++)
          if(imol = molecule[j]) {

/* uncomment maybe? */
//        if (molmap) imol = molmap[imol-idlo];
//        else imol--;
            domain->unmap(x[j],image[j],unwrap);
            if (rmass) massone = rmass[j];
            massone = mass[type[j]];		

            x[i][0] += unwrap[0] * massone;
            x[i][1] += unwrap[1] * massone;
            x[i][2] += unwrap[2] * massone;
            masstotal += massone; 
  
    }
     x[i][0] /= masstotal;
     x[i][1] /= masstotal;
     x[i][2] /= masstotal;      
}

//  MPI_Allreduce(&com[0][0],&comall[0][0],3*nmolecules,
//                MPI_DOUBLE,MPI_SUM,world);
//  for (int i = 0; i < nmolecules; i++) {
//    comall[i][0] /= masstotal[i];
//    comall[i][1] /= masstotal[i];
 //   comall[i][2] /= masstotal[i];

}

/* ---------------------------------------------------------------------- */

void FixDistForce::min_post_force(int vflag)
{
  post_force(vflag);
}
