// This file was auto-generated

%include "fp_common.i"

%{
#include "vlidort_interface_masters.h"
%}

%import "vlidort_interface_types.i"

%fp_shared_ptr(FullPhysics::VBrdf_Linsup_Masters);
%fp_shared_ptr(FullPhysics::VBrdf_Sup_Masters);
%fp_shared_ptr(FullPhysics::VLidort_Inputs);
%fp_shared_ptr(FullPhysics::VLidort_Masters);
%fp_shared_ptr(FullPhysics::VLidort_L_Inputs);
%fp_shared_ptr(FullPhysics::VLidort_Lcs_Masters);
%fp_shared_ptr(FullPhysics::VLidort_Lps_Masters);
%fp_shared_ptr(FullPhysics::VLidort_Vbrdf_Sup_Accessories);

namespace FullPhysics {



class VBrdf_Linsup_Masters {

public:
  VBrdf_Linsup_Masters();
  virtual ~VBrdf_Linsup_Masters();
  std::string print_to_string() const;

  %python_attribute(vbrdf_sup_in, VBrdf_Sup_Inputs&)
  %python_attribute(vbrdf_linsup_in, VBrdf_Linsup_Inputs&)
  %python_attribute(vbrdf_sup_inputstatus, VBrdf_Input_Exception_Handling&)
  %python_attribute(vbrdf_sup_out, VBrdf_Sup_Outputs&)
  %python_attribute(vbrdf_linsup_out, VBrdf_Linsup_Outputs&)
  %python_attribute(vbrdf_sup_outputstatus, VBrdf_Output_Exception_Handling&)
  
  void read_config(const std::string& filnam_in);
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in);
};


class VBrdf_Sup_Masters {

public:
  VBrdf_Sup_Masters();
  virtual ~VBrdf_Sup_Masters();
  std::string print_to_string() const;

  %python_attribute(vbrdf_sup_in, VBrdf_Sup_Inputs&)
  %python_attribute(vbrdf_sup_inputstatus, VBrdf_Input_Exception_Handling&)
  %python_attribute(vbrdf_sup_out, VBrdf_Sup_Outputs&)
  %python_attribute(vbrdf_sup_outputstatus, VBrdf_Output_Exception_Handling&)
  
  void read_config(const std::string& filnam_in);
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in);
};


class VLidort_Inputs {

public:
  VLidort_Inputs();
  virtual ~VLidort_Inputs();
  std::string print_to_string() const;

  %python_attribute(vlidort_sup, VLidort_Sup_Inout&)
  %python_attribute(vlidort_fixin, VLidort_Fixed_Inputs&)
  %python_attribute(vlidort_modin, VLidort_Modified_Inputs&)
  %python_attribute(vlidort_inputstatus, VLidort_Input_Exception_Handling&)
  
  void brdf_sup_init();
  void read_config(const std::string& filnam_in);
  void sleave_sup_init();
  void ss_sup_init();
  void sup_init();
};


class VLidort_Masters {

public:
  VLidort_Masters();
  virtual ~VLidort_Masters();
  std::string print_to_string() const;

  %python_attribute(vlidort_fixin, VLidort_Fixed_Inputs&)
  %python_attribute(vlidort_modin, VLidort_Modified_Inputs&)
  %python_attribute(vlidort_sup, VLidort_Sup_Inout&)
  %python_attribute(vlidort_out, VLidort_Outputs&)
  
  void run(const bool& do_debug_input_in);
};


class VLidort_L_Inputs {

public:
  VLidort_L_Inputs();
  virtual ~VLidort_L_Inputs();
  std::string print_to_string() const;

  %python_attribute(vlidort_linsup, VLidort_Linsup_Inout&)
  %python_attribute(vlidort_fixin, VLidort_Fixed_Inputs&)
  %python_attribute(vlidort_modin, VLidort_Modified_Inputs&)
  %python_attribute(vlidort_linfixin, VLidort_Fixed_Lininputs&)
  %python_attribute(vlidort_linmodin, VLidort_Modified_Lininputs&)
  %python_attribute(vlidort_inputstatus, VLidort_Input_Exception_Handling&)
  
  void brdf_linsup_init();
  void read_config(const std::string& filnam_in);
  void linsup_init();
  void sleave_linsup_init();
  void ss_linsup_init();
};


class VLidort_Lcs_Masters {

public:
  VLidort_Lcs_Masters();
  virtual ~VLidort_Lcs_Masters();
  std::string print_to_string() const;

  %python_attribute(vlidort_fixin, VLidort_Fixed_Inputs&)
  %python_attribute(vlidort_modin, VLidort_Modified_Inputs&)
  %python_attribute(vlidort_sup, VLidort_Sup_Inout&)
  %python_attribute(vlidort_out, VLidort_Outputs&)
  %python_attribute(vlidort_linfixin, VLidort_Fixed_Lininputs&)
  %python_attribute(vlidort_linmodin, VLidort_Modified_Lininputs&)
  %python_attribute(vlidort_linsup, VLidort_Linsup_Inout&)
  %python_attribute(vlidort_linout, VLidort_Linoutputs&)
  
  void run(const bool& do_debug_input_in);
};


class VLidort_Lps_Masters {

public:
  VLidort_Lps_Masters();
  virtual ~VLidort_Lps_Masters();
  std::string print_to_string() const;

  %python_attribute(vlidort_fixin, VLidort_Fixed_Inputs&)
  %python_attribute(vlidort_modin, VLidort_Modified_Inputs&)
  %python_attribute(vlidort_sup, VLidort_Sup_Inout&)
  %python_attribute(vlidort_out, VLidort_Outputs&)
  %python_attribute(vlidort_linfixin, VLidort_Fixed_Lininputs&)
  %python_attribute(vlidort_linmodin, VLidort_Modified_Lininputs&)
  %python_attribute(vlidort_linsup, VLidort_Linsup_Inout&)
  %python_attribute(vlidort_linout, VLidort_Linoutputs&)
  
  void run(const bool& do_debug_input_in);
};


class VLidort_Vbrdf_Sup_Accessories {

public:
  VLidort_Vbrdf_Sup_Accessories();
  virtual ~VLidort_Vbrdf_Sup_Accessories();
  std::string print_to_string() const;

  %python_attribute(vbrdf_sup_out, VBrdf_Sup_Outputs&)
  %python_attribute(vlidort_fixin, VLidort_Fixed_Inputs&)
  %python_attribute(vlidort_modin, VLidort_Modified_Inputs&)
  %python_attribute(vlidort_sup, VLidort_Sup_Inout&)
  %python_attribute(vbrdf_sup_in, VBrdf_Sup_Inputs&)
  %python_attribute(vlidort_vbrdfcheck_status, VLidort_Exception_Handling&)
  
  void set_vbrdf_inputs();
  void vbrdf_input_check();
  void vbrdf_input_check_error(const std::string& errorfile_in);
};

}