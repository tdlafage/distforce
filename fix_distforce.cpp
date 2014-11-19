// fix fix_id group_id(this should be the virtual particles) distforce

#include "fix_distforce.h"
#include "error.h"
#include "atom.h"
#include "memory.h"
#include "group.h"
#include "domain.h"
#include "mpi.h"
#include "update.h"
#include "compute.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixDistForce::FixDistForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

   int nmax=atom->nmax;
   tagint idlo;
   tagint idhi;

// Error flags. First one makes sure the correct number of flags. The second
// makes sure that there molecule numbers in the atom type.

   if (narg != 3) error->all(FLERR,"Illegal fix distforce command");

   if (atom->molecular == 0)
    error->all(FLERR,"Fix distforce requires full or molecular atom style");

   nevery=1;

  nmolecules = 250;
  size_array_rows = nmolecules;

// Memory allocation

   memory->create(masstotal,nmolecules,"distforce/molecule:masstotal");
   memory->create(massproc,nmolecules,"distforce/molecule:massproc");
   memory->create(comall,nmolecules,3,"distforce/molecule:comall");
   memory->create(comhold,nmolecules,3,"distforce/molecule:comhold");
   meomry->create(forceproc,nmolecules,3,"distforce/molecule:forceproc");
   memory->create(forceall,nmolecules,3,"distforce/molecule:forceall");

  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  tagint imol;
  double massone;

// Molecular mass calculator

  for (int i = 0; i < nmolecules; i++) massproc[i] = 0.0;

  for (int i = 0; i < nlocal; i++){
	imol = molecule[i];
	if (rmass) massone = rmass[i];
	else massone = mass[type[i]];
	if(!(mask[i] & groupbit)){
		massproc[imol-1] += massone;
	}
  }

  MPI_Allreduce(&massproc[0],&masstotal[0],nmolecules,MPI_DOUBLE,MPI_SUM,world);

}

/* ---------------------------------------------------------------------- */

FixDistForce::~FixDistForce()
{
  memory->destroy(masstotal);
  memory->destroy(massproc);
  memory->destroy(comall);
  memory->destroy(comhold);
  memory->destroy(forceproc);
  memory->destroy(forceall);
}

/* ---------------------------------------------------------------------- */
using namespace FixConst;
int FixDistForce::setmask()
{
  // set with a bitmask how and when apply the force from EDM 
  int mask = 0;
  mask |=  FixConst::POST_FORCE;
  mask |=  FixConst::MIN_POST_FORCE;
  mask |=  FixConst::POST_INTEGRATE;
  mask |=  FixConst::MIN_PRE_EXCHANGE; 
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

void FixDistForce::post_force(int vflag)
{
  tagint imol;
  double massone;

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

// Force distributor

  for(int i = 0; i < nmolecules; i++){
	forceproc[i][0] = 0.0;
	forceproc[i][1] = 0.0;
	forceproc[i][2] = 0.0;
  }

  for(int i = 0; i < nlocal; i++){
	if(mask[i] & groubit){
		imol = molecule[i];
		forceproc[imol-1][0] = f[i][0];
		forceproc[imol-1][1] = f[i][1];
		forceproc[imol-1][2] = f[i][2];

		f[i][0] = 0.0;
		f[i][1] = 0.0;
		f[i][2] = 0.0;
	}
  }

  MPI_Allreduce(&forceproc[0][0],&forceall[0][0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);

  for(int i = 0; i < nlocal; i++){
	if(!(mask[i] & groupbit){
		imol = molecule[i];
		if (rmass) massone = rmass[i];
		else massone = mass[type[i]];
		f[i][0] += massone*forceall[imol-1][0]/masstotal[imol-1];
		f[i][1] += massone*forceall[imol-1][1]/masstotal[imol-1];
		f[i][2] += massone*forceall[imol-1][2]/masstotal[imol-1];
	}
  }
}

  

/* ---------------------------------------------------------------------- */

void FixDistForce::post_integrate()
 {
  tagint imol;
  double massone;
  double unwrap[3];
  int n;
  int k;


  double *lo = domain->boxlo;
  double *hi = domain->boxhi;
  double **x = atom->x;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;
  tagint *tag = atom->tag;

// COM pseudo integrator

  for(int i = 0; i < nmolecules; i++){
	comhold[i][0] = 0.0;
	comhold[i][1] = 0.0;
	comhold[i][2] = 0.0;
  }

  for(int i = 0; i < nlocal; i++){
	imol = molecule[i];
	if (rmass) massone = rmass[i];
	else massone = mass[type[i]];
	if(!(mask[i] & groupbit)){
		domain->unmap(x[i],image[i],unwrap);
		comhold[imol-1][0] += unwrap[0] * massone / masstotal[imol-1];
		comhold[imol-1][1] += unwrap[1] * massone / masstotal[imol-1];
		comhold[imol-1][2] += unwrap[2] * massone / masstotal[imol-1];
	}
  }

  MPI_Allreduce(&comhold[0][0],&comall[0][0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);


  for(int i = 0; i < nlocal; i++){
	if(mask[i] & groupbit) {
		imol = molecule[i];
		imol--;
		x[i][0] = comall[imol][0];
		x[i][1] = comall[imol][1];
		x[i][2] = comall[imol][2];

//		image flag calculator

		for(int k = 0;k < 3;k++){
			if(x[i][k] > hi[k]){
				n = int((x[i][k]+lo[k])/(hi[k]-lo[k]));
				image[i] = n;
			}
			if(x[i][k] < lo[k]){
				n = int((x[i][k]-hi[k])/(hi[k]-lo[k]));
				image[i] = n;
			}
		}
	}
  }
}


/* ---------------------------------------------------------------------- */

void FixDistForce::min_post_force(int vflag)
{
}

/* ---------------------------------------------------------------------- */

void FixDistForce::min_pre_exchange()
{
}

