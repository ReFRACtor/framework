#ifndef MAX_A_POSTERIORI_SQRT_CONSTRAINT_H
#define MAX_A_POSTERIORI_SQRT_CONSTRAINT_H
#include "max_a_posteriori.h"
#include "model_measure_standard.h"
#include "state_mapping_linear.h"
#include "observer.h"
#include <boost/make_shared.hpp>

namespace FullPhysics {
/******************************************************************
  This is a variation of MaxAPosterioriStandard, where we supply
  a "sqrt_constraint" to use instead of the a_priori matrix.

  In addition, we take a StateMapping that maps between our "retrieval
  vector" and "full state vector". We can support most of the
  functionality by just using StateMapping in the various pieces that
  make up the StateVector, however muses actually needs access to
  jacobian on the full forward model grid, K_x. So we need to track
  this in this class, the forward model generates K_x and we then
  handle calculating K_z. For a discussion of this, see section
  III.A.1 of "Tropospheric Emission Spectrometer: Retrieval Method and
  Error Analysis" (IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING,
  VOL. 44, NO. 5, MAY 2006).

  Our original MaxAPosterioriStandard uses a Cholesky decomposition,
  while py-retrieve uses a different decomposition based on SVD
  (there are numerous decomposition of B = A A^T).
  
  The two are in some sense equivalent - the norm of the residual
  for the apriori augmentation is the same. However the specific
  values of the residual are different - basically these are rotations
  of each other.

  I'm not sure that it matters, but there may be calculation in
  py-retrieve that depends on the specific order of the residual
  relative to its sqrt_constraint matrix. So we have this variation
  that uses a supplied sqrt_constraint matrix rather than taking an 
  a priori matrix. Note that the sqrt_constraint we take is the
  transpose of what py-retrieve has - for some reason it transposes
  this after calculating it.

  Note by the way that despite the name, sqrt_constraint *isn't*
  actually the sqrt of the constraint matrix. The

  constraint_matrix = a_priori_cov^-1

  And the matrix is decomposed so

  sqrt_constraint * sqrt_constraint^T = constraint_matrix

  The decomposition by something like scipy.linalg.sqrtm is also
  a common way to handle the a priori augmentation, but this has

  r * r = constraint_matrix

  if r = sqrtm(constraint_matrix).

  sqrt_constraint is *not* equal to r, py-retrieve doesn't actually
  do a sqrtm.

  Instead the sqrt_constraint is like a cholesky decomposition of the
  constraint_matrix, it just isn't confined to being a lower
  triangular matrix. I don't believe the decomposition is unique
  either (although I'm not sure about that) - but since we take the
  decomposition as an input there isn't any ambiguity in this class.

  Note if you use this class then the calls to a_priori_cov_chol_inv()
  and a_priori_cov_chol() actually returns sqrt_constraint and
  sqrt_constraint^-1. I briefly considered overriding these functions
  to raise an Exception, since these aren't actually lower triangular
  matrices you get from a Cholesky decomposition. But I decide against
  that, for most uses I could think of you probably don't actually
  care that this is a different decomposition. But we can reconsider
  that if needed in the future - we really do have function names here
  that lie.
*******************************************************************/
class MaxAPosterioriSqrtConstraint : 
    virtual public MaxAPosteriori, virtual public ModelMeasureStandard,
    public Observable<MaxAPosterioriSqrtConstraint> {
public:
  MaxAPosterioriSqrtConstraint(const boost::shared_ptr<ForwardModel>& fm,
          const boost::shared_ptr<Observation>& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint,
          const boost::shared_ptr<StateMapping>& in_map = boost::make_shared<StateMappingLinear>());

  MaxAPosterioriSqrtConstraint(const std::vector<boost::shared_ptr<ForwardModel> >& fm,
          const std::vector<boost::shared_ptr<Observation> >& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint,
          const boost::shared_ptr<StateMapping>& in_map =
          boost::make_shared<StateMappingLinear>());
  
  virtual ~MaxAPosterioriSqrtConstraint() {}

  virtual void add_observer(Observer<MaxAPosterioriSqrtConstraint>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<MaxAPosterioriSqrtConstraint>& Obs) 
  { remove_observer_do(Obs, *this);}
  
//-----------------------------------------------------------------------
/// Jacobian on full state grid (K_x in nomenclature of TES paper)
//-----------------------------------------------------------------------
  
  virtual blitz::Array<double, 2> jacobian_fm()
  { jacobian_eval(); return K_x.copy(); }

//-----------------------------------------------------------------------
/// Measurment jacobian on full state grid (K_x in nomenclature of TES paper)
//-----------------------------------------------------------------------
  
  virtual blitz::Array<double, 2> measurement_jacobian_fm()
  { measurement_jacobian_eval(); return msrmnt_jacobian_x.copy(); }
  
//-----------------------------------------------------------------------
/// Jacobian of model - measurement on full state grid.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> model_measure_diff_jacobian_fm();
  
//-----------------------------------------------------------------------
/// StateMapping to go to and from
//-----------------------------------------------------------------------

  const boost::shared_ptr<StateMapping>& mapping() const
  { return mapping_;}

  virtual void parameters(const blitz::Array<double, 1>& x);
  virtual blitz::Array<double, 1> parameters() const
  { return ModelMeasureStandard::parameters(); }
  
  virtual int expected_parameter_size() const { return Xa.rows(); }

  void get_state(bool& msrmnt_is_const_v,
		 blitz::Array<double, 1>& M_v,
		 blitz::Array<double, 2>& K_v,
		 blitz::Array<double, 1>& msrmnt_v,
		 blitz::Array<double, 2>& msrmnt_jacobian_v,
		 blitz::Array<double, 2>& K_x_v,
		 blitz::Array<double, 2>& msrmnt_jacobian_x_v) const;
  
  void set_state(const bool& msrmnt_is_const_v,
		 const blitz::Array<double, 1>& M_v,
		 const blitz::Array<double, 2>& K_v,
		 const blitz::Array<double, 1>& msrmnt_v,
		 const blitz::Array<double, 2>& msrmnt_jacobian_v,
		 const blitz::Array<double, 2>& K_x_v,
		 const blitz::Array<double, 2>& msrmnt_jacobian_x_v);
		 
		 
//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "MaxAPosterioriSqrtConstraint"; }

protected:
  virtual void radiance_from_fm(bool skip_check=false);
  virtual void measurement_eval();
private:
  blitz::Array<double, 2> K_x;
  blitz::Array<double, 2> msrmnt_jacobian_x;
  boost::shared_ptr<StateMapping> mapping_;
  MaxAPosterioriSqrtConstraint() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MaxAPosterioriSqrtConstraint);
FP_EXPORT_OBSERVER_KEY(MaxAPosterioriSqrtConstraint);
#endif
