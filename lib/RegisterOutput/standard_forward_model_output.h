#ifndef OCO_FORWARD_MODEL_OUTPUT_H
#define OCO_FORWARD_MODEL_OUTPUT_H
#include "register_output_base.h"
#include "standard_forward_model.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the StandardForwardModel class that 
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just 
  part of the StandardForwardModel class.
*******************************************************************/
class StandardForwardModelOutput : public RegisterOutputBase {
public:
    StandardForwardModelOutput(const boost::shared_ptr<StandardForwardModel>&
                       Fm)
        : fm(Fm) {}
    virtual ~StandardForwardModelOutput() {}
    virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
    boost::shared_ptr<StandardForwardModel> fm;
};
}
#endif
