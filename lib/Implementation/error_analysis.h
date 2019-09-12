#ifndef ERROR_ANALYSIS_H
#define ERROR_ANALYSIS_H
#include "connor_solver.h"
#include "max_a_posteriori.h"
#include "atmosphere_standard.h"
#include "forward_model.h"
#include "observation.h"
#include "fe_disable_exception.h"

namespace FullPhysics {
/****************************************************************//**
  This calculates a variety of values to help with the error analysis
  of a Level 2 Full Physics Run.

  We currently support both a ConnorSolver or the more general
  MaxAPosteriori. The error analysis is almost identical, we just
  get parameters from one or the other source.

  Note that the current implementation of this class repeatedly
  calculates certain values (e.g., hmat is calculated each time it is
  used). We could cache values if needed, with code handling the
  clearing of the cache when absorber or solver changes. But right now
  the error analysis is only done a hand full of times (in a normal
  run, just at the end. With iteration output, at each
  iteration). Performance is perfectly acceptable, even with duplicate
  calculations. We can revisit this if performance ever becomes an
  issue.

  \todo The various component calculation in section 3.6.5 and 3.6.6
  in the ATB assume that the state vector contains the mixing ratio of
  the CO2 on levels. What needs to change in these calculations if it
  doesn't? (e.g., a scale, a shape)
*******************************************************************/

class ErrorAnalysis : public Printable<ErrorAnalysis> {
public:
    ErrorAnalysis(const boost::shared_ptr<ConnorSolver>& Solver,
                  const boost::shared_ptr<RtAtmosphere>& Atm,
                  const boost::shared_ptr<ForwardModel>& Fm,
                  const boost::shared_ptr<Observation>& inst_meas);
    ErrorAnalysis(const boost::shared_ptr<MaxAPosteriori>& Max_a_posteriori,
                  const boost::shared_ptr<RtAtmosphere>& Atm,
                  const boost::shared_ptr<ForwardModel>& Fm,
                  const boost::shared_ptr<Observation>& inst_meas);
    virtual ~ErrorAnalysis() {}

    int num_channels() const;

    blitz::Array<double, 1> residual() const;
    blitz::Array<double, 1> modeled_radiance() const;

    double signal_level(int band) const;
    double noise_level(int band) const;

    blitz::Array<double, 2> aposteriori_covariance() const;
    blitz::Array<double, 2> apriori_covariance() const;

    blitz::Array<double, 2> averaging_kernel() const;

    double residual_sum_sq(int Band) const;
    double residual_mean_sq(int Band) const;

    double reduced_chisq(int Band) const;

    double relative_residual_mean_sq(int Band) const;

    double chisq_measure_norm(const blitz::Array<double, 1>& Residual,
                              const blitz::Array<double, 1>& Residual_cov_diag) const;

    double degrees_of_freedom() const;

    void print(std::ostream& Os) const
    {
        Os << "ErrorAnalysis";
    }
private:

    // Only one of solver or max_a_posteriori will be nonnull.
    boost::shared_ptr<ConnorSolver> solver;
    boost::shared_ptr<MaxAPosteriori> max_a_posteriori;
    boost::shared_ptr<AtmosphereStandard> atm;
    boost::shared_ptr<ForwardModel> fm;
    boost::shared_ptr<Observation> meas;

    // Used in a lot of places, so define once here.
    blitz::firstIndex i1;
    blitz::secondIndex i2;
    blitz::thirdIndex i3;
    blitz::fourthIndex i4;
};
}
#endif
