// This file was auto-generated

%include "fp_common.i"

%{
#include "spurr_interface_masters.h"
%}

%import "spurr_interface_types.i"

%fp_shared_ptr(FullPhysics::Spurr_Brdf_Lin_Sup_Masters_Base);
%fp_shared_ptr(FullPhysics::Spurr_Brdf_Sup_Masters_Base);
%fp_shared_ptr(FullPhysics::Spurr_Inputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Masters_Base);
%fp_shared_ptr(FullPhysics::Spurr_L_Inputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Lcs_Masters_Base);
%fp_shared_ptr(FullPhysics::Spurr_Lps_Masters_Base);
%fp_shared_ptr(FullPhysics::Spurr_Brdf_Sup_Accessories_Base);

namespace FullPhysics {


class Spurr_Brdf_Lin_Sup_Masters_Base {

public:
  std::string print_to_string() const;

  %python_attribute_abstract(brdf_sup_inputs_base, Spurr_Brdf_Sup_Inputs_Base&)
  %python_attribute_abstract(brdf_linsup_inputs_base, Spurr_Brdf_Lin_Sup_Inputs_Base&)
  %python_attribute_abstract(brdf_input_exception_handling_base, Spurr_Brdf_Input_Exception_Handling_Base&)
  %python_attribute_abstract(brdf_sup_outputs_base, Spurr_Brdf_Sup_Outputs_Base&)
  %python_attribute_abstract(brdf_linsup_outputs_base, Spurr_Brdf_Lin_Sup_Outputs_Base&)
  %python_attribute_abstract(brdf_output_exception_handling_base, Spurr_Brdf_Output_Exception_Handling_Base&)
  
  void read_config(const std::string& filnam_in) = 0;
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) = 0;
};

class Spurr_Brdf_Sup_Masters_Base {

public:
  std::string print_to_string() const;

  %python_attribute_abstract(brdf_sup_inputs_base, Spurr_Brdf_Sup_Inputs_Base&)
  %python_attribute_abstract(brdf_input_exception_handling_base, Spurr_Brdf_Input_Exception_Handling_Base&)
  %python_attribute_abstract(brdf_sup_outputs_base, Spurr_Brdf_Sup_Outputs_Base&)
  %python_attribute_abstract(brdf_output_exception_handling_base, Spurr_Brdf_Output_Exception_Handling_Base&)
  
  void read_config(const std::string& filnam_in) = 0;
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) = 0;
};

class Spurr_Inputs_Base {

public:
  std::string print_to_string() const;

  %python_attribute_abstract(sup_inout_base, Spurr_Sup_Inout_Base&)
  %python_attribute_abstract(fixed_inputs_base, Spurr_Fixed_Inputs_Base&)
  %python_attribute_abstract(modified_inputs_base, Spurr_Modified_Inputs_Base&)
  %python_attribute_abstract(input_exception_handling_base, Spurr_Input_Exception_Handling_Base&)
  
  void read_config(const std::string& filnam_in) = 0;
};

class Spurr_Masters_Base {

public:
  std::string print_to_string() const;

  %python_attribute_abstract(fixed_inputs_base, Spurr_Fixed_Inputs_Base&)
  %python_attribute_abstract(modified_inputs_base, Spurr_Modified_Inputs_Base&)
  %python_attribute_abstract(sup_inout_base, Spurr_Sup_Inout_Base&)
  %python_attribute_abstract(outputs_base, Spurr_Outputs_Base&)
  
  void run(const bool& do_debug_input_in) = 0;
};

class Spurr_L_Inputs_Base {

public:
  std::string print_to_string() const;

  %python_attribute_abstract(linsup_inout_base, Spurr_Lin_Sup_Inout_Base&)
  %python_attribute_abstract(fixed_inputs_base, Spurr_Fixed_Inputs_Base&)
  %python_attribute_abstract(modified_inputs_base, Spurr_Modified_Inputs_Base&)
  %python_attribute_abstract(fixed_lininputs_base, Spurr_Fixed_Lininputs_Base&)
  %python_attribute_abstract(modified_lininputs_base, Spurr_Modified_Lininputs_Base&)
  %python_attribute_abstract(input_exception_handling_base, Spurr_Input_Exception_Handling_Base&)
  
  void read_config(const std::string& filnam_in) = 0;
};

class Spurr_Lcs_Masters_Base {

public:
  std::string print_to_string() const;

  %python_attribute_abstract(fixed_inputs_base, Spurr_Fixed_Inputs_Base&)
  %python_attribute_abstract(modified_inputs_base, Spurr_Modified_Inputs_Base&)
  %python_attribute_abstract(sup_inout_base, Spurr_Sup_Inout_Base&)
  %python_attribute_abstract(outputs_base, Spurr_Outputs_Base&)
  %python_attribute_abstract(fixed_lininputs_base, Spurr_Fixed_Lininputs_Base&)
  %python_attribute_abstract(modified_lininputs_base, Spurr_Modified_Lininputs_Base&)
  %python_attribute_abstract(linsup_inout_base, Spurr_Lin_Sup_Inout_Base&)
  %python_attribute_abstract(linoutputs_base, Spurr_Linoutputs_Base&)
  
  void run(const bool& do_debug_input_in) = 0;
};

class Spurr_Lps_Masters_Base {

public:
  std::string print_to_string() const;

  %python_attribute_abstract(fixed_inputs_base, Spurr_Fixed_Inputs_Base&)
  %python_attribute_abstract(modified_inputs_base, Spurr_Modified_Inputs_Base&)
  %python_attribute_abstract(sup_inout_base, Spurr_Sup_Inout_Base&)
  %python_attribute_abstract(outputs_base, Spurr_Outputs_Base&)
  %python_attribute_abstract(fixed_lininputs_base, Spurr_Fixed_Lininputs_Base&)
  %python_attribute_abstract(modified_lininputs_base, Spurr_Modified_Lininputs_Base&)
  %python_attribute_abstract(linsup_inout_base, Spurr_Lin_Sup_Inout_Base&)
  %python_attribute_abstract(linoutputs_base, Spurr_Linoutputs_Base&)
  
  void run(const bool& do_debug_input_in) = 0;
};

class Spurr_Brdf_Sup_Accessories_Base {

public:
  std::string print_to_string() const;

  %python_attribute_abstract(brdf_sup_inputs_base, Spurr_Brdf_Sup_Inputs_Base&)
  %python_attribute_abstract(fixed_inputs_base, Spurr_Fixed_Inputs_Base&)
  %python_attribute_abstract(modified_inputs_base, Spurr_Modified_Inputs_Base&)
  %python_attribute_abstract(exception_handling_base, Spurr_Exception_Handling_Base&)
  %python_attribute_abstract(brdf_sup_outputs_base, Spurr_Brdf_Sup_Outputs_Base&)
  %python_attribute_abstract(sup_inout_base, Spurr_Sup_Inout_Base&)
  
};

}