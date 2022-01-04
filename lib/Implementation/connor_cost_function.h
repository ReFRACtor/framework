#ifndef CONNOR_COST_FUNCTION
#define CONNOR_COST_FUNCTION
#include "cost_function.h"
#include "state_vector.h"
#include "forward_model.h"
#include "observation.h"
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
class ConnorCostFunction : public CostFunction {
public:
  ConnorCostFunction(
          const boost::shared_ptr<StateVector>& Sv, 
          const boost::shared_ptr<ForwardModel>& fm, 
          const boost::shared_ptr<Observation>& inst_meas)
    : statev(Sv), forward_model(fm), meas(inst_meas)
  {
  }
  virtual ~ConnorCostFunction() {}
  virtual void cost_function(const blitz::Array<double, 1>& X,
                        blitz::Array<double, 1>& Residual,
                        blitz::Array<double, 1>& Se,
                        blitz::Array<double, 2>& Jacobian) const
  {
    using namespace blitz;
    statev->update_state(X);
    SpectralRange rad_meas = meas->radiance_all().spectral_range();
    Spectrum rad_spec = forward_model->radiance_all();

    SpectralRange rad_mod = rad_spec.spectral_range().convert(rad_meas.units());
    Se.resize(rad_meas.data().rows());
    Residual.resize(rad_meas.data().rows());

    if(rad_meas.uncertainty().rows() == 0)
      throw Exception("Radiance uncertainty is empty in ConnorCostFunction");

    Se = sqr(rad_meas.uncertainty());
    Residual = rad_mod.data() - rad_meas.data();
    Jacobian.reference(rad_mod.data_ad().jacobian());
  }
  const boost::shared_ptr<Observation> observation() const { return meas; }
  virtual void print(std::ostream& Os) const
  {
    Os << "ConnorCostFunction";
  }
private:
  boost::shared_ptr<StateVector> statev;
  boost::shared_ptr<ForwardModel> forward_model;
  boost::shared_ptr<Observation> meas;

  ConnorCostFunction() = default;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(ConnorCostFunction);

#endif
