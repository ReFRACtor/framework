#include <max_a_posteriori_standard.h>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(MaxAPosterioriStandard, MaxAPosteriori)
.def(luabind::constructor< const boost::shared_ptr<ForwardModel>&,
                           const boost::shared_ptr<Observation>&, 
                           const boost::shared_ptr<StateVector>&,
                           const Array<double, 1>,
                           const Array<double, 2> >())
REGISTER_LUA_END()
#endif



MaxAPosterioriStandard::MaxAPosterioriStandard(const boost::shared_ptr<ForwardModel>& forward_model,
        const boost::shared_ptr<Observation>& observation, 
        const boost::shared_ptr<StateVector>& state_vector,
        const Array<double, 1> a_priori_params,
        const Array<double, 2> a_priori_cov)
  : ModelMeasure(
        obs->radiance_all().spectral_range().data(),
        Array<double, 1>(sqr(obs->radiance_all().spectral_range().uncertainty()))),
    MaxAPosteriori(a_priori_params, a_priori_cov),
    ModelMeasureStandard(forward_model, observation, state_vector) 
{
  if(Xa.rows() != sv->observer_claimed_size())
    throw Exception("A priori state vector size and state vector size expected by the model are not equal. :( ");
}


//  TEMPORARY
//
// Should go away after we end support for 
// fixed pressure level grid. 
// Not implemented efficiently.
#include <linear_algebra.h>
void MaxAPosterioriStandard::vanishing_params_update()
{
  ModelMeasureStandard::vanishing_params_update();
  Array<bool, 1> used(sv->used_flag());
  if(used.rows() != Sa_chol.rows())
    throw Exception("Size of a-priori cov. matrix and the number of the elements of used-flag inconsistent! :( ");
  if(!all(used)) {
//    throw Exception("Handling vanishing parameters is not supported yet! :( ");
    Sa_chol.reference(cholesky_decomposition(Sa));
    for(int i=0; i<used.rows(); i++)
      if(!used(i))
        Sa_chol(i,Range::all()) = 0.0;
    Sa_chol_inv.reference(generalized_inverse(Sa_chol,1e-20));
    // Theoretically the selected columns of Sa_chol_inv should
    // be zero; however, they may have very small nonzero values.
    for(int i=0; i<used.rows(); i++)
      if(!used(i))
        Sa_chol_inv(Range::all(),i) = 0.0;
  }
}
