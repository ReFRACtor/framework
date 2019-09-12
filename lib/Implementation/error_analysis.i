%include "fp_common.i"
%{
#include "error_analysis.h"
#include "sub_state_vector_array.h"
%}
%base_import(generic_object)
%import "connor_solver.i"
%import "max_a_posteriori.i"
%import "atmosphere_standard.i"
%import "forward_model.i"
%import "observation.i"
%fp_shared_ptr(FullPhysics::ErrorAnalysis);

namespace FullPhysics {
class ErrorAnalysis : public GenericObject {
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

    %python_attribute(degrees_of_freedom, double)
};
}
