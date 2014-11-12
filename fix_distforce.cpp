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

   int nmax=atom->nmax;

   if (narg != 3) error->all(FLERR,"Illegal fix distforce command");

   if (atom->molecular == 0)
    error->all(FLERR,"Fix distforce requires full or molecular atom style");

   nevery=1;

   memory->create(masstotal,nmax,"distforce/molecule:masstotal");
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
//  mask |= PRE_FORCE;
//  mask |= MIN_PRE_FORCE;
  mask |=  FixConst::POST_FORCE;
  mask |=  FixConst::MIN_POST_FORCE;
  mask |=  FixConst::END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDistForce::init()
{

   int nmax=atom->nmax;

  for(int i = 0; i < nmax; i++){
	masstotal[i] = -1.0;
  }

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

//void FixDistForce::pre_force(int vflag)
//{
//}
/* ---------------------------------------------------------------------- */

void FixDistForce::post_force(int vflag)
{
  int imol;
  double massone;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  for (int i = 0; i < nlocal; i++){
     if (mask[i] & groupbit) {
	imol = molecule[i];
	for(int k = 0; k < nlocal; k++){
	if(imol == molecule[k]){
	if(i != k){
		if(masstotal[imol] <  0.0){
		summolmass(imol);
		}
//		if (rmass) massone = rmass[k];
		massone = mass[type[k]];
		f[k][0] += massone*f[i][0]/masstotal[imol];
		f[k][1] += massone*f[i][1]/masstotal[imol];
		f[k][2] += massone*f[i][2]/masstotal[imol];
	}
	}
	}
   }
   }
}

/* ---------------------------------------------------------------------- */

void FixDistForce::end_of_step()
 {
  int imol;
  double massone;
  double unwrap[3];

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


  for (int i = 0; i < (nlocal); i++){
	if (mask[i] & groupbit) {
		x[i][0] = x[i][1] = x[i][2] = 0.0;
		imol = molecule[i];
//    if (molmap) imol = molmap[imol-idlo];
//    else imol--;
		for(int j = 0; j < (nlocal); j++){
			if(imol == molecule[j]) {
				if(i != j){

/* uncomment maybe? */
//        if (molmap) imol = molmap[imol-idlo];
//        else imol--;
					domain->unmap(x[j],image[j],unwrap);
//					if (rmass) massone = rmass[j];
					massone = mass[type[j]];		

					x[i][0] += unwrap[0] * massone;
					x[i][1] += unwrap[1] * massone;
					x[i][2] += unwrap[2] * massone;
				}
			}
		}
		if(masstotal[imol] < 0.0){
			summolmass(imol);
		}
//		printf("%lf\n",masstotal[imol]);
		x[i][0] /= masstotal[imol];
		x[i][1] /= masstotal[imol];
		x[i][2] /= masstotal[imol];
	}
    }
      
}

//  MPI_Allreduce(&com[0][0],&comall[0][0],3*nmolecules,
//                MPI_DOUBLE,MPI_SUM,world);
//  for (int i = 0; i < nmolecules; i++) {
//    comall[i][0] /= masstotal[i];
//    comall[i][1] /= masstotal[i];
//    comall[i][2] /= masstotal[i];

/* ---------------------------------------------------------------------- */

void FixDistForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixDistForce::summolmass(int num)
{
	double massone;

	double *mass = atom->mass;
	double *rmass = atom->rmass;
	int nlocal = atom->nlocal;
	int nghost = atom->nghost;
	tagint *molecule = atom->molecule;
	int *type = atom->type;

	int imol;

//	printf("summolmass is being called");
//	printf("nlocal:%d nghosts:%d\n",nlocal,nghost);
	masstotal[num] = 0.0;

	for(int i = 0; i < (nlocal); i++){
		imol = molecule[i];
		if(imol == num){
			massone = mass[type[i]];
			masstotal[num] += massone;
		}
	}
//	printf("Molecular mass:%lf\n",masstotal[num]);		
}
