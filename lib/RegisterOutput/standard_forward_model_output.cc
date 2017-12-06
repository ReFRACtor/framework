#include "standard_forward_model_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> oco_fm_out_create
(const boost::shared_ptr<ForwardModel>& fm)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new StandardForwardModelOutput
     (boost::dynamic_pointer_cast<StandardForwardModel>(fm)));
}
REGISTER_LUA_DERIVED_CLASS(StandardForwardModelOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &oco_fm_out_create)
]
REGISTER_LUA_END()
#endif

// Adapter class to take things from StandardForwardModel into format needed
// for output
class StandardForwardModelCalc {
public:
  StandardForwardModelCalc(const boost::shared_ptr<StandardForwardModel>& Fm)
    : fm(Fm) {}

  Array<int, 1> samples() const 
  {
    int num_total_samples = 0;
    for(int i = 0; i < fm->num_channels(); ++i) {
        num_total_samples += fm->spectral_domain(i).data().rows();
    }
    Array<int, 1> samples(num_total_samples);
    for(int i = 0; i < fm->num_channels(); ++i) {
      boost::optional<blitz::Range> pr = fm->stacked_pixel_range(i);
      if(pr)
        samples(*pr) = fm->spectral_grid()->low_resolution_grid(i).sample_index();
    }
    return samples;
  }

private:
    boost::shared_ptr<StandardForwardModel> fm;
};

// See base class for description

void StandardForwardModelOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    boost::shared_ptr<StandardForwardModelCalc> fmc(new StandardForwardModelCalc(fm));
    out->register_data_source("/SpectralParameters/sample_indexes",
                              &StandardForwardModelCalc::samples, fmc);
}
