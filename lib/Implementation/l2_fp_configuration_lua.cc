#include "l2_fp_configuration_lua.h"
#include "hdf_file_generating.h"
#include "output_hdf.h"
#include "output_hdf_iteration.h"
#include "state_vector.h"
using namespace FullPhysics;

// See base class for description.
void L2FpConfigurationLua::output(boost::shared_ptr<Output>& Regular_output,
				  boost::shared_ptr<Output>& Error_output) const
{
  boost::shared_ptr<HdfFileGenerating> output_g(new HdfFileGenerating(output_name_));
  int num_aer_part = ls->globals()["number_aerosol"].value<int>();
  int npres = ls->globals()["number_pressure_level"].value<int>();
  int nband = ls->globals()["number_band"].value<int>();
  int sv_claimed_size = ls->globals()["state_vector"].value_ptr<StateVector>()->observer_claimed_size();
  Regular_output.reset
    (new OutputHdf(output_g, npres, 
		   sv_claimed_size, 
		   num_aer_part + 1, nband));
  if(ls->globals()["iteration_output"].value<bool>()) {
    if(solver()) {
      boost::shared_ptr<OutputHdfIteration> outit(new OutputHdfIteration(output_g));
      out_iteration = outit;
      solver()->add_observer(*outit);
    } else
      throw Exception("Iteration_output can only be used with a Connor solver, not the newer IterativeSolver that Edwin developed");
  }
  Error_output.reset
    (new OutputHdf(output_name_ + ".error", npres, 
		   sv_claimed_size, 
		   num_aer_part + 1, nband));

  /* Removed output, previously done via RegisterOutput family here */
}
