#ifndef SPECTRAL_PARAMETERS_OUTPUT_H
#define SPECTRAL_PARAMETERS_OUTPUT_H

#include "register_output_base.h"
#include "forward_model.h"
#include "observation.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the ForwardModel class that should be
  written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the ForwardModel class.
*******************************************************************/
class SpectralParametersOutput : public RegisterOutputBase {
public:
    SpectralParametersOutput(const boost::shared_ptr<ForwardModel>& Fm, const boost::shared_ptr<Observation>& inst_meas)
        : fm(Fm), meas(inst_meas) {}
    virtual ~SpectralParametersOutput() {}
    virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
    boost::shared_ptr<ForwardModel> fm;
    boost::shared_ptr<Observation> meas;
};
}
#endif
