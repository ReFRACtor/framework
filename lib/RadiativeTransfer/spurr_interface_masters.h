#ifndef SPURR_INTERFACE_MASTERS_H
#define SPURR_INTERFACE_MASTERS_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"
#include "spurr_interface_types.h"


/* This file was auto-generated */

namespace FullPhysics {

class Spurr_Brdf_Lin_Sup_Masters_Base : public virtual Printable<Spurr_Brdf_Lin_Sup_Masters_Base> {

public:

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  virtual Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() = 0;
  virtual const Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() const = 0;
  
  virtual Spurr_Brdf_Lin_Sup_Inputs_Base& brdf_linsup_inputs_base() = 0;
  virtual const Spurr_Brdf_Lin_Sup_Inputs_Base& brdf_linsup_inputs_base() const = 0;
  
  virtual Spurr_Brdf_Input_Exception_Handling_Base& brdf_input_exception_handling_base() = 0;
  virtual const Spurr_Brdf_Input_Exception_Handling_Base& brdf_input_exception_handling_base() const = 0;
  
  virtual Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() = 0;
  virtual const Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() const = 0;
  
  virtual Spurr_Brdf_Lin_Sup_Outputs_Base& brdf_linsup_outputs_base() = 0;
  virtual const Spurr_Brdf_Lin_Sup_Outputs_Base& brdf_linsup_outputs_base() const = 0;
  
  virtual Spurr_Brdf_Output_Exception_Handling_Base& brdf_output_exception_handling_base() = 0;
  virtual const Spurr_Brdf_Output_Exception_Handling_Base& brdf_output_exception_handling_base() const = 0;
  
  
  virtual void read_config(const std::string& filnam_in) = 0;
  virtual void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;

};

class Spurr_Brdf_Sup_Masters_Base : public virtual Printable<Spurr_Brdf_Sup_Masters_Base> {

public:

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  virtual Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() = 0;
  virtual const Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() const = 0;
  
  virtual Spurr_Brdf_Input_Exception_Handling_Base& brdf_input_exception_handling_base() = 0;
  virtual const Spurr_Brdf_Input_Exception_Handling_Base& brdf_input_exception_handling_base() const = 0;
  
  virtual Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() = 0;
  virtual const Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() const = 0;
  
  virtual Spurr_Brdf_Output_Exception_Handling_Base& brdf_output_exception_handling_base() = 0;
  virtual const Spurr_Brdf_Output_Exception_Handling_Base& brdf_output_exception_handling_base() const = 0;
  
  
  virtual void read_config(const std::string& filnam_in) = 0;
  virtual void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;

};

class Spurr_Inputs_Base : public virtual Printable<Spurr_Inputs_Base> {

public:

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  virtual Spurr_Sup_Inout_Base& sup_inout_base() = 0;
  virtual const Spurr_Sup_Inout_Base& sup_inout_base() const = 0;
  
  virtual Spurr_Fixed_Inputs_Base& fixed_inputs_base() = 0;
  virtual const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const = 0;
  
  virtual Spurr_Modified_Inputs_Base& modified_inputs_base() = 0;
  virtual const Spurr_Modified_Inputs_Base& modified_inputs_base() const = 0;
  
  virtual Spurr_Input_Exception_Handling_Base& input_exception_handling_base() = 0;
  virtual const Spurr_Input_Exception_Handling_Base& input_exception_handling_base() const = 0;
  
  
  virtual void read_config(const std::string& filnam_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;

};

class Spurr_Masters_Base : public virtual Printable<Spurr_Masters_Base> {

public:

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  virtual Spurr_Fixed_Inputs_Base& fixed_inputs_base() = 0;
  virtual const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const = 0;
  
  virtual Spurr_Modified_Inputs_Base& modified_inputs_base() = 0;
  virtual const Spurr_Modified_Inputs_Base& modified_inputs_base() const = 0;
  
  virtual Spurr_Sup_Inout_Base& sup_inout_base() = 0;
  virtual const Spurr_Sup_Inout_Base& sup_inout_base() const = 0;
  
  virtual Spurr_Outputs_Base& outputs_base() = 0;
  virtual const Spurr_Outputs_Base& outputs_base() const = 0;
  
  
  virtual void run(const bool& do_debug_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;

};

class Spurr_L_Inputs_Base : public virtual Printable<Spurr_L_Inputs_Base> {

public:

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  virtual Spurr_Lin_Sup_Inout_Base& linsup_inout_base() = 0;
  virtual const Spurr_Lin_Sup_Inout_Base& linsup_inout_base() const = 0;
  
  virtual Spurr_Fixed_Inputs_Base& fixed_inputs_base() = 0;
  virtual const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const = 0;
  
  virtual Spurr_Modified_Inputs_Base& modified_inputs_base() = 0;
  virtual const Spurr_Modified_Inputs_Base& modified_inputs_base() const = 0;
  
  virtual Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() = 0;
  virtual const Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() const = 0;
  
  virtual Spurr_Modified_Lininputs_Base& modified_lininputs_base() = 0;
  virtual const Spurr_Modified_Lininputs_Base& modified_lininputs_base() const = 0;
  
  virtual Spurr_Input_Exception_Handling_Base& input_exception_handling_base() = 0;
  virtual const Spurr_Input_Exception_Handling_Base& input_exception_handling_base() const = 0;
  
  
  virtual void read_config(const std::string& filnam_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;

};

class Spurr_Lcs_Masters_Base : public virtual Printable<Spurr_Lcs_Masters_Base> {

public:

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  virtual Spurr_Fixed_Inputs_Base& fixed_inputs_base() = 0;
  virtual const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const = 0;
  
  virtual Spurr_Modified_Inputs_Base& modified_inputs_base() = 0;
  virtual const Spurr_Modified_Inputs_Base& modified_inputs_base() const = 0;
  
  virtual Spurr_Sup_Inout_Base& sup_inout_base() = 0;
  virtual const Spurr_Sup_Inout_Base& sup_inout_base() const = 0;
  
  virtual Spurr_Outputs_Base& outputs_base() = 0;
  virtual const Spurr_Outputs_Base& outputs_base() const = 0;
  
  virtual Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() = 0;
  virtual const Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() const = 0;
  
  virtual Spurr_Modified_Lininputs_Base& modified_lininputs_base() = 0;
  virtual const Spurr_Modified_Lininputs_Base& modified_lininputs_base() const = 0;
  
  virtual Spurr_Lin_Sup_Inout_Base& linsup_inout_base() = 0;
  virtual const Spurr_Lin_Sup_Inout_Base& linsup_inout_base() const = 0;
  
  virtual Spurr_Linoutputs_Base& linoutputs_base() = 0;
  virtual const Spurr_Linoutputs_Base& linoutputs_base() const = 0;
  
  
  virtual void run(const bool& do_debug_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;

};

class Spurr_Lps_Masters_Base : public virtual Printable<Spurr_Lps_Masters_Base> {

public:

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  virtual Spurr_Fixed_Inputs_Base& fixed_inputs_base() = 0;
  virtual const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const = 0;
  
  virtual Spurr_Modified_Inputs_Base& modified_inputs_base() = 0;
  virtual const Spurr_Modified_Inputs_Base& modified_inputs_base() const = 0;
  
  virtual Spurr_Sup_Inout_Base& sup_inout_base() = 0;
  virtual const Spurr_Sup_Inout_Base& sup_inout_base() const = 0;
  
  virtual Spurr_Outputs_Base& outputs_base() = 0;
  virtual const Spurr_Outputs_Base& outputs_base() const = 0;
  
  virtual Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() = 0;
  virtual const Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() const = 0;
  
  virtual Spurr_Modified_Lininputs_Base& modified_lininputs_base() = 0;
  virtual const Spurr_Modified_Lininputs_Base& modified_lininputs_base() const = 0;
  
  virtual Spurr_Lin_Sup_Inout_Base& linsup_inout_base() = 0;
  virtual const Spurr_Lin_Sup_Inout_Base& linsup_inout_base() const = 0;
  
  virtual Spurr_Linoutputs_Base& linoutputs_base() = 0;
  virtual const Spurr_Linoutputs_Base& linoutputs_base() const = 0;
  
  
  virtual void run(const bool& do_debug_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;

};

class Spurr_Brdf_Sup_Accessories_Base : public virtual Printable<Spurr_Brdf_Sup_Accessories_Base> {

public:

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  virtual Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() = 0;
  virtual const Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() const = 0;
  
  virtual Spurr_Fixed_Inputs_Base& fixed_inputs_base() = 0;
  virtual const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const = 0;
  
  virtual Spurr_Modified_Inputs_Base& modified_inputs_base() = 0;
  virtual const Spurr_Modified_Inputs_Base& modified_inputs_base() const = 0;
  
  virtual Spurr_Exception_Handling_Base& exception_handling_base() = 0;
  virtual const Spurr_Exception_Handling_Base& exception_handling_base() const = 0;
  
  virtual Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() = 0;
  virtual const Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() const = 0;
  
  virtual Spurr_Sup_Inout_Base& sup_inout_base() = 0;
  virtual const Spurr_Sup_Inout_Base& sup_inout_base() const = 0;
  
  
  
  virtual void print(std::ostream &output_stream) const = 0;

};

}
#endif