

#include "fix_distforce.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixDistForce::FixDistForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

FixDistForce::~FixDistForce()
{
  
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
}

/* ---------------------------------------------------------------------- */

void FixDistForce::min_post_force(int vflag)
{
  post_force(vflag);
}
