#include "state_mapping.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateMapping::serialize(Archive& ar, const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(StateMapping);

  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
}

FP_IMPLEMENT(StateMapping);
#endif

blitz::Array<double, 2> StateMapping::jacobian_retrieval
  (const blitz::Array<double, 1>& retrieval_vector,
   const blitz::Array<double, 2>& jacobian_mapped) const
{
  // In the nomenclature of the TES paper section
  // III.A.1 of "Tropospheric Emission Spectrometer: Retrieval Method and
  // Error Analysis" (IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING,
  // VOL. 44, NO. 5, MAY 2006)
  // The retrieval_vector is "z", the mapped full state vector is "x"
  // jacobian_mapped is K_x (df/dx). We get dz/dx for our current
  // state, and calculate K_z (df/dz) as df/dx * dx/dz
  blitz::firstIndex i1; blitz::secondIndex i2; blitz::thirdIndex i3;
  blitz::Array<double, 2> dz_dz(retrieval_vector.rows(), retrieval_vector.rows());
  dz_dz = 0;
  for(int i = 0; i < dz_dz.rows(); ++i)
    dz_dz(i,i) = 1;
  blitz::Array<double, 2> dx_dz = mapped_state(ArrayAd<double, 1>(retrieval_vector, dz_dz)).jacobian();
  if(jacobian_mapped.cols() != dx_dz.rows()) {
    Exception e("jacobian_mapped cols needs to match the mapped vector rows\n");
    e << "jacobian_mapped shape: " << jacobian_mapped.rows() << " x "
      << jacobian_mapped.cols() << "\n"
      << "dx_dz shape: " << dx_dz.rows() << " x "
      << dx_dz.cols() << "\n";
    throw e;
  }
  blitz::Array<double, 2> res(jacobian_mapped.rows(), dx_dz.cols());
  res = blitz::sum(jacobian_mapped(i1, i3) * dx_dz(i3, i2), i3);
  return res;
}
