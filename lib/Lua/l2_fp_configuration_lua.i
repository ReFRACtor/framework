// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "l2_fp_configuration_lua.h"
#include "state_vector.h"
%}
%base_import(l2_fp_configuration)
%import "lua_state.i"
%import "state_vector.i"
%import "output.i"
%fp_shared_ptr(FullPhysics::L2FpConfigurationLua);

namespace FullPhysics {
class LuaState;
class L2FpConfigurationLua : public L2FpConfiguration {
public:
  L2FpConfigurationLua(const std::string& Fname, 
		       const std::string& Out_file = "out.h5");
  L2FpConfigurationLua(const boost::shared_ptr<LuaState>& Ls, 
		       const std::string& Out_file = "out.h5");
  L2FpConfigurationLua(int Argc, char** Argv);
  %python_attribute_nonconst(lua_state, LuaState)
  %python_attribute_with_set(output_name, std::string);
  %python_attribute(state_vector, boost::shared_ptr<StateVector>)
  virtual void output(boost::shared_ptr<Output>& OUTPUT,
		      boost::shared_ptr<Output>& OUTPUT) const;
};
}
