#include "convergence_check.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ConvergenceCheck::serialize(Archive & ar,
                                 const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(ConvergenceCheck);

  // Dummy placeholder, just so we can have derived classes call
  // serialization of this. We use to have derived classes "know"
  // that the base class doesn't have anything. But seems better to
  // *always* have base classes do something, so we can add stuff in
  // the future w/o breaking a bunch of code.
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
}

template<class Archive>
void FitStatistic::serialize(Archive & ar,
                             const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(FitStatistic);
  ar & FP_NVP(fit_succeeded) & FP_NVP(outcome)
    & FP_NVP(number_iteration) & FP_NVP(number_divergent)
    & FP_NVP(d_sigma_sq) & FP_NVP(d_sigma_sq_scaled)
    & FP_NVP(chisq_apriori) & FP_NVP(chisq_measured)
    & FP_NVP(chisq_apriori_fc) & FP_NVP(chisq_measured_fc);
}

FP_IMPLEMENT(ConvergenceCheck);
FP_IMPLEMENT(FitStatistic)
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(ConvergenceCheck)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void FitStatistic::print(std::ostream& Os) const
{
  Os << "FitStatistic\n"
     << "  number_iteration:  " << number_iteration << "\n"
     << "  number_divergent:  " << number_divergent << "\n"
     << "  outcome:           " << (int) outcome << "\n"
     << "  fit_succeeded:     " << fit_succeeded << "\n"
     << "  d_sigma_sq:        " << d_sigma_sq << "\n"
     << "  d_sigma_sq_scaled: " << d_sigma_sq_scaled << "\n"
     << "  chisq_apriori:     " << chisq_apriori << "\n"
     << "  chisq_measured:    " << chisq_measured << "\n"
     << "  chisq_apriori_fc:  " << chisq_apriori_fc << "\n"
     << "  chisq_measured_fc: " << chisq_measured_fc << "\n"
     << "  gamma2:            " << gamma2() << "\n"
     << "  gamma2_fc:         " << gamma2_fc() << "\n";
}

