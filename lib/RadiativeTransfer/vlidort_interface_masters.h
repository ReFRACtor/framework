#ifndef VLIDORT_INTERFACE_MASTERS_H
#define VLIDORT_INTERFACE_MASTERS_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"
#include "vlidort_interface_types.h"
#include "spurr_interface_masters.h"


/* This file was auto-generated */

namespace FullPhysics {

//-----------------------------------------------------------------------
// Links to module: "vbrdf_linsup_masters_m" in file: "vbrdf_lin_sup_masters.f90"
//-----------------------------------------------------------------------

extern "C" {
  void vbrdf_linsup_masters_m_read_wrap(const char* filename, const int* filename_len, void** vbrdf_sup_in_in, void** vbrdf_linsup_in_in, void** vbrdf_sup_inputstatus_in, void** vbrdf_sup_out_in, void** vbrdf_linsup_out_in, void** vbrdf_sup_outputstatus_in);
  void vbrdf_linsup_masters_m_write_wrap(const char* filename, const int* filename_len, void** vbrdf_sup_in_in, void** vbrdf_linsup_in_in, void** vbrdf_sup_inputstatus_in, void** vbrdf_sup_out_in, void** vbrdf_linsup_out_in, void** vbrdf_sup_outputstatus_in);
  void vbrdf_linsup_masters_m_vbrdf_lin_inputmaster_wrap(const int* filnam_in_len, const char* filnam_in, void** vbrdf_sup_in_in, void** vbrdf_linsup_in_in, void** vbrdf_sup_inputstatus_in);
  void vbrdf_linsup_masters_m_vbrdf_lin_mainmaster_wrap(const bool* do_debug_restoration_in, const int* nmoments_input_in, void** vbrdf_sup_in_in, void** vbrdf_linsup_in_in, void** vbrdf_sup_out_in, void** vbrdf_linsup_out_in, void** vbrdf_sup_outputstatus_in);
}

class VBrdf_Linsup_Masters : public virtual Spurr_Brdf_Lin_Sup_Masters_Base {

public:
  VBrdf_Linsup_Masters() 
  { 
    
    // Initialize type pointers
    vbrdf_sup_in_.reset( new VBrdf_Sup_Inputs() );
    vbrdf_linsup_in_.reset( new VBrdf_Linsup_Inputs() );
    vbrdf_sup_inputstatus_.reset( new VBrdf_Input_Exception_Handling() );
    vbrdf_sup_out_.reset( new VBrdf_Sup_Outputs() );
    vbrdf_linsup_out_.reset( new VBrdf_Linsup_Outputs() );
    vbrdf_sup_outputstatus_.reset( new VBrdf_Output_Exception_Handling() );
    
  }

  virtual ~VBrdf_Linsup_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  VBrdf_Sup_Inputs& vbrdf_sup_in() {
    return *vbrdf_sup_in_;
  }

  const VBrdf_Sup_Inputs& vbrdf_sup_in() const {
    return *vbrdf_sup_in_;
  }

  Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() {
    return *vbrdf_sup_in_;
  }

  const Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() const {
    return *vbrdf_sup_in_;
  }

  boost::shared_ptr<VBrdf_Sup_Inputs>& vbrdf_sup_in_ptr() {
    return vbrdf_sup_in_;
  }

  
  VBrdf_Linsup_Inputs& vbrdf_linsup_in() {
    return *vbrdf_linsup_in_;
  }

  const VBrdf_Linsup_Inputs& vbrdf_linsup_in() const {
    return *vbrdf_linsup_in_;
  }

  Spurr_Brdf_Lin_Sup_Inputs_Base& brdf_linsup_inputs_base() {
    return *vbrdf_linsup_in_;
  }

  const Spurr_Brdf_Lin_Sup_Inputs_Base& brdf_linsup_inputs_base() const {
    return *vbrdf_linsup_in_;
  }

  boost::shared_ptr<VBrdf_Linsup_Inputs>& vbrdf_linsup_in_ptr() {
    return vbrdf_linsup_in_;
  }

  
  VBrdf_Input_Exception_Handling& vbrdf_sup_inputstatus() {
    return *vbrdf_sup_inputstatus_;
  }

  const VBrdf_Input_Exception_Handling& vbrdf_sup_inputstatus() const {
    return *vbrdf_sup_inputstatus_;
  }

  Spurr_Brdf_Input_Exception_Handling_Base& brdf_input_exception_handling_base() {
    return *vbrdf_sup_inputstatus_;
  }

  const Spurr_Brdf_Input_Exception_Handling_Base& brdf_input_exception_handling_base() const {
    return *vbrdf_sup_inputstatus_;
  }

  boost::shared_ptr<VBrdf_Input_Exception_Handling>& vbrdf_sup_inputstatus_ptr() {
    return vbrdf_sup_inputstatus_;
  }

  
  VBrdf_Sup_Outputs& vbrdf_sup_out() {
    return *vbrdf_sup_out_;
  }

  const VBrdf_Sup_Outputs& vbrdf_sup_out() const {
    return *vbrdf_sup_out_;
  }

  Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() {
    return *vbrdf_sup_out_;
  }

  const Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() const {
    return *vbrdf_sup_out_;
  }

  boost::shared_ptr<VBrdf_Sup_Outputs>& vbrdf_sup_out_ptr() {
    return vbrdf_sup_out_;
  }

  
  VBrdf_Linsup_Outputs& vbrdf_linsup_out() {
    return *vbrdf_linsup_out_;
  }

  const VBrdf_Linsup_Outputs& vbrdf_linsup_out() const {
    return *vbrdf_linsup_out_;
  }

  Spurr_Brdf_Lin_Sup_Outputs_Base& brdf_linsup_outputs_base() {
    return *vbrdf_linsup_out_;
  }

  const Spurr_Brdf_Lin_Sup_Outputs_Base& brdf_linsup_outputs_base() const {
    return *vbrdf_linsup_out_;
  }

  boost::shared_ptr<VBrdf_Linsup_Outputs>& vbrdf_linsup_out_ptr() {
    return vbrdf_linsup_out_;
  }

  
  VBrdf_Output_Exception_Handling& vbrdf_sup_outputstatus() {
    return *vbrdf_sup_outputstatus_;
  }

  const VBrdf_Output_Exception_Handling& vbrdf_sup_outputstatus() const {
    return *vbrdf_sup_outputstatus_;
  }

  Spurr_Brdf_Output_Exception_Handling_Base& brdf_output_exception_handling_base() {
    return *vbrdf_sup_outputstatus_;
  }

  const Spurr_Brdf_Output_Exception_Handling_Base& brdf_output_exception_handling_base() const {
    return *vbrdf_sup_outputstatus_;
  }

  boost::shared_ptr<VBrdf_Output_Exception_Handling>& vbrdf_sup_outputstatus_ptr() {
    return vbrdf_sup_outputstatus_;
  }

  
  void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vbrdf_linsup_in_lcl = vbrdf_linsup_in_->fortran_type_ptr();
    void* vbrdf_sup_inputstatus_lcl = vbrdf_sup_inputstatus_->fortran_type_ptr();
    
    vbrdf_linsup_masters_m_vbrdf_lin_inputmaster_wrap(&filnam_in_len, filnam_lcl, &vbrdf_sup_in_lcl, &vbrdf_linsup_in_lcl, &vbrdf_sup_inputstatus_lcl);
    

  }

  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) {
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vbrdf_linsup_in_lcl = vbrdf_linsup_in_->fortran_type_ptr();
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vbrdf_linsup_out_lcl = vbrdf_linsup_out_->fortran_type_ptr();
    void* vbrdf_sup_outputstatus_lcl = vbrdf_sup_outputstatus_->fortran_type_ptr();
    
    vbrdf_linsup_masters_m_vbrdf_lin_mainmaster_wrap(&do_debug_restoration_in, &nmoments_input_in, &vbrdf_sup_in_lcl, &vbrdf_linsup_in_lcl, &vbrdf_sup_out_lcl, &vbrdf_linsup_out_lcl, &vbrdf_sup_outputstatus_lcl);
    

  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vbrdf_linsup_in_lcl = vbrdf_linsup_in_->fortran_type_ptr();
    void* vbrdf_sup_inputstatus_lcl = vbrdf_sup_inputstatus_->fortran_type_ptr();
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vbrdf_linsup_out_lcl = vbrdf_linsup_out_->fortran_type_ptr();
    void* vbrdf_sup_outputstatus_lcl = vbrdf_sup_outputstatus_->fortran_type_ptr();

    vbrdf_linsup_masters_m_read_wrap(filename_lcl, &filename_in_len, &vbrdf_sup_in_lcl, &vbrdf_linsup_in_lcl, &vbrdf_sup_inputstatus_lcl, &vbrdf_sup_out_lcl, &vbrdf_linsup_out_lcl, &vbrdf_sup_outputstatus_lcl);
    
  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) const {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vbrdf_linsup_in_lcl = vbrdf_linsup_in_->fortran_type_ptr();
    void* vbrdf_sup_inputstatus_lcl = vbrdf_sup_inputstatus_->fortran_type_ptr();
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vbrdf_linsup_out_lcl = vbrdf_linsup_out_->fortran_type_ptr();
    void* vbrdf_sup_outputstatus_lcl = vbrdf_sup_outputstatus_->fortran_type_ptr();

    vbrdf_linsup_masters_m_write_wrap(filename_lcl, &filename_in_len, &vbrdf_sup_in_lcl, &vbrdf_linsup_in_lcl, &vbrdf_sup_inputstatus_lcl, &vbrdf_sup_out_lcl, &vbrdf_linsup_out_lcl, &vbrdf_sup_outputstatus_lcl);
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "VBrdf_Linsup_Masters:" << std::endl
      << "          vbrdf_sup_in: " << vbrdf_sup_in()  << std::endl
      << "       vbrdf_linsup_in: " << vbrdf_linsup_in()  << std::endl
      << " vbrdf_sup_inputstatus: " << vbrdf_sup_inputstatus()  << std::endl
      << "         vbrdf_sup_out: " << vbrdf_sup_out()  << std::endl
      << "      vbrdf_linsup_out: " << vbrdf_linsup_out()  << std::endl
      << "vbrdf_sup_outputstatus: " << vbrdf_sup_outputstatus()  << std::endl;

  }

private:
  boost::shared_ptr<VBrdf_Sup_Inputs> vbrdf_sup_in_;
  boost::shared_ptr<VBrdf_Linsup_Inputs> vbrdf_linsup_in_;
  boost::shared_ptr<VBrdf_Input_Exception_Handling> vbrdf_sup_inputstatus_;
  boost::shared_ptr<VBrdf_Sup_Outputs> vbrdf_sup_out_;
  boost::shared_ptr<VBrdf_Linsup_Outputs> vbrdf_linsup_out_;
  boost::shared_ptr<VBrdf_Output_Exception_Handling> vbrdf_sup_outputstatus_;

  // Serialization support
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "vbrdf_sup_masters_m" in file: "vbrdf_sup_masters.f90"
//-----------------------------------------------------------------------

extern "C" {
  void vbrdf_sup_masters_m_read_wrap(const char* filename, const int* filename_len, void** vbrdf_sup_in_in, void** vbrdf_sup_inputstatus_in, void** vbrdf_sup_out_in, void** vbrdf_sup_outputstatus_in);
  void vbrdf_sup_masters_m_write_wrap(const char* filename, const int* filename_len, void** vbrdf_sup_in_in, void** vbrdf_sup_inputstatus_in, void** vbrdf_sup_out_in, void** vbrdf_sup_outputstatus_in);
  void vbrdf_sup_masters_m_vbrdf_inputmaster_wrap(const int* filnam_in_len, const char* filnam_in, void** vbrdf_sup_in_in, void** vbrdf_sup_inputstatus_in);
  void vbrdf_sup_masters_m_vbrdf_mainmaster_wrap(const bool* do_debug_restoration_in, const int* nmoments_input_in, void** vbrdf_sup_in_in, void** vbrdf_sup_out_in, void** vbrdf_sup_outputstatus_in);
}

class VBrdf_Sup_Masters : public virtual Spurr_Brdf_Sup_Masters_Base {

public:
  VBrdf_Sup_Masters() 
  { 
    
    // Initialize type pointers
    vbrdf_sup_in_.reset( new VBrdf_Sup_Inputs() );
    vbrdf_sup_inputstatus_.reset( new VBrdf_Input_Exception_Handling() );
    vbrdf_sup_out_.reset( new VBrdf_Sup_Outputs() );
    vbrdf_sup_outputstatus_.reset( new VBrdf_Output_Exception_Handling() );
    
  }

  virtual ~VBrdf_Sup_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  VBrdf_Sup_Inputs& vbrdf_sup_in() {
    return *vbrdf_sup_in_;
  }

  const VBrdf_Sup_Inputs& vbrdf_sup_in() const {
    return *vbrdf_sup_in_;
  }

  Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() {
    return *vbrdf_sup_in_;
  }

  const Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() const {
    return *vbrdf_sup_in_;
  }

  boost::shared_ptr<VBrdf_Sup_Inputs>& vbrdf_sup_in_ptr() {
    return vbrdf_sup_in_;
  }

  
  VBrdf_Input_Exception_Handling& vbrdf_sup_inputstatus() {
    return *vbrdf_sup_inputstatus_;
  }

  const VBrdf_Input_Exception_Handling& vbrdf_sup_inputstatus() const {
    return *vbrdf_sup_inputstatus_;
  }

  Spurr_Brdf_Input_Exception_Handling_Base& brdf_input_exception_handling_base() {
    return *vbrdf_sup_inputstatus_;
  }

  const Spurr_Brdf_Input_Exception_Handling_Base& brdf_input_exception_handling_base() const {
    return *vbrdf_sup_inputstatus_;
  }

  boost::shared_ptr<VBrdf_Input_Exception_Handling>& vbrdf_sup_inputstatus_ptr() {
    return vbrdf_sup_inputstatus_;
  }

  
  VBrdf_Sup_Outputs& vbrdf_sup_out() {
    return *vbrdf_sup_out_;
  }

  const VBrdf_Sup_Outputs& vbrdf_sup_out() const {
    return *vbrdf_sup_out_;
  }

  Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() {
    return *vbrdf_sup_out_;
  }

  const Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() const {
    return *vbrdf_sup_out_;
  }

  boost::shared_ptr<VBrdf_Sup_Outputs>& vbrdf_sup_out_ptr() {
    return vbrdf_sup_out_;
  }

  
  VBrdf_Output_Exception_Handling& vbrdf_sup_outputstatus() {
    return *vbrdf_sup_outputstatus_;
  }

  const VBrdf_Output_Exception_Handling& vbrdf_sup_outputstatus() const {
    return *vbrdf_sup_outputstatus_;
  }

  Spurr_Brdf_Output_Exception_Handling_Base& brdf_output_exception_handling_base() {
    return *vbrdf_sup_outputstatus_;
  }

  const Spurr_Brdf_Output_Exception_Handling_Base& brdf_output_exception_handling_base() const {
    return *vbrdf_sup_outputstatus_;
  }

  boost::shared_ptr<VBrdf_Output_Exception_Handling>& vbrdf_sup_outputstatus_ptr() {
    return vbrdf_sup_outputstatus_;
  }

  
  void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vbrdf_sup_inputstatus_lcl = vbrdf_sup_inputstatus_->fortran_type_ptr();
    
    vbrdf_sup_masters_m_vbrdf_inputmaster_wrap(&filnam_in_len, filnam_lcl, &vbrdf_sup_in_lcl, &vbrdf_sup_inputstatus_lcl);
    

  }

  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) {
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vbrdf_sup_outputstatus_lcl = vbrdf_sup_outputstatus_->fortran_type_ptr();
    
    vbrdf_sup_masters_m_vbrdf_mainmaster_wrap(&do_debug_restoration_in, &nmoments_input_in, &vbrdf_sup_in_lcl, &vbrdf_sup_out_lcl, &vbrdf_sup_outputstatus_lcl);
    

  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vbrdf_sup_inputstatus_lcl = vbrdf_sup_inputstatus_->fortran_type_ptr();
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vbrdf_sup_outputstatus_lcl = vbrdf_sup_outputstatus_->fortran_type_ptr();

    vbrdf_sup_masters_m_read_wrap(filename_lcl, &filename_in_len, &vbrdf_sup_in_lcl, &vbrdf_sup_inputstatus_lcl, &vbrdf_sup_out_lcl, &vbrdf_sup_outputstatus_lcl);
    
  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) const {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vbrdf_sup_inputstatus_lcl = vbrdf_sup_inputstatus_->fortran_type_ptr();
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vbrdf_sup_outputstatus_lcl = vbrdf_sup_outputstatus_->fortran_type_ptr();

    vbrdf_sup_masters_m_write_wrap(filename_lcl, &filename_in_len, &vbrdf_sup_in_lcl, &vbrdf_sup_inputstatus_lcl, &vbrdf_sup_out_lcl, &vbrdf_sup_outputstatus_lcl);
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "VBrdf_Sup_Masters:" << std::endl
      << "          vbrdf_sup_in: " << vbrdf_sup_in()  << std::endl
      << " vbrdf_sup_inputstatus: " << vbrdf_sup_inputstatus()  << std::endl
      << "         vbrdf_sup_out: " << vbrdf_sup_out()  << std::endl
      << "vbrdf_sup_outputstatus: " << vbrdf_sup_outputstatus()  << std::endl;

  }

private:
  boost::shared_ptr<VBrdf_Sup_Inputs> vbrdf_sup_in_;
  boost::shared_ptr<VBrdf_Input_Exception_Handling> vbrdf_sup_inputstatus_;
  boost::shared_ptr<VBrdf_Sup_Outputs> vbrdf_sup_out_;
  boost::shared_ptr<VBrdf_Output_Exception_Handling> vbrdf_sup_outputstatus_;

  // Serialization support
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "vlidort_inputs_m" in file: "vlidort_inputs.f90"
//-----------------------------------------------------------------------

extern "C" {
  void v_inputs_m_read_wrap(const char* filename, const int* filename_len, void** vlidort_sup_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_inputstatus_in);
  void v_inputs_m_write_wrap(const char* filename, const int* filename_len, void** vlidort_sup_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_inputstatus_in);
  void v_inputs_m_v_brdf_sup_init_wrap(void** vlidort_sup_in);
  void v_inputs_m_v_input_master_wrap(const int* filnam_in_len, const char* filnam_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_inputstatus_in);
  void v_inputs_m_v_sleave_sup_init_wrap(void** vlidort_sup_in);
  void v_inputs_m_v_ss_sup_init_wrap(void** vlidort_sup_in);
  void v_inputs_m_v_sup_init_wrap(void** vlidort_sup_in);
}

class VLidort_Inputs : public virtual Spurr_Inputs_Base {

public:
  VLidort_Inputs() 
  { 
    
    // Initialize type pointers
    vlidort_sup_.reset( new VLidort_Sup_Inout() );
    vlidort_fixin_.reset( new VLidort_Fixed_Inputs() );
    vlidort_modin_.reset( new VLidort_Modified_Inputs() );
    vlidort_inputstatus_.reset( new VLidort_Input_Exception_Handling() );
    
  }

  virtual ~VLidort_Inputs() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  VLidort_Sup_Inout& vlidort_sup() {
    return *vlidort_sup_;
  }

  const VLidort_Sup_Inout& vlidort_sup() const {
    return *vlidort_sup_;
  }

  Spurr_Sup_Inout_Base& sup_inout_base() {
    return *vlidort_sup_;
  }

  const Spurr_Sup_Inout_Base& sup_inout_base() const {
    return *vlidort_sup_;
  }

  boost::shared_ptr<VLidort_Sup_Inout>& vlidort_sup_ptr() {
    return vlidort_sup_;
  }

  void vlidort_sup(VLidort_Sup_Inout& vlidort_sup_in) {
    void* src_ptr = vlidort_sup_in.fortran_type_ptr();
    void* dst_ptr = vlidort_sup_->fortran_type_ptr();
    vlidort_sup_inout_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Fixed_Inputs& vlidort_fixin() {
    return *vlidort_fixin_;
  }

  const VLidort_Fixed_Inputs& vlidort_fixin() const {
    return *vlidort_fixin_;
  }

  Spurr_Fixed_Inputs_Base& fixed_inputs_base() {
    return *vlidort_fixin_;
  }

  const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const {
    return *vlidort_fixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Inputs>& vlidort_fixin_ptr() {
    return vlidort_fixin_;
  }

  
  VLidort_Modified_Inputs& vlidort_modin() {
    return *vlidort_modin_;
  }

  const VLidort_Modified_Inputs& vlidort_modin() const {
    return *vlidort_modin_;
  }

  Spurr_Modified_Inputs_Base& modified_inputs_base() {
    return *vlidort_modin_;
  }

  const Spurr_Modified_Inputs_Base& modified_inputs_base() const {
    return *vlidort_modin_;
  }

  boost::shared_ptr<VLidort_Modified_Inputs>& vlidort_modin_ptr() {
    return vlidort_modin_;
  }

  
  VLidort_Input_Exception_Handling& vlidort_inputstatus() {
    return *vlidort_inputstatus_;
  }

  const VLidort_Input_Exception_Handling& vlidort_inputstatus() const {
    return *vlidort_inputstatus_;
  }

  Spurr_Input_Exception_Handling_Base& input_exception_handling_base() {
    return *vlidort_inputstatus_;
  }

  const Spurr_Input_Exception_Handling_Base& input_exception_handling_base() const {
    return *vlidort_inputstatus_;
  }

  boost::shared_ptr<VLidort_Input_Exception_Handling>& vlidort_inputstatus_ptr() {
    return vlidort_inputstatus_;
  }

  
  void brdf_sup_init() {
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    
    v_inputs_m_v_brdf_sup_init_wrap(&vlidort_sup_lcl);
    

  }

  void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_inputstatus_lcl = vlidort_inputstatus_->fortran_type_ptr();
    
    v_inputs_m_v_input_master_wrap(&filnam_in_len, filnam_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_inputstatus_lcl);
    

  }

  void sleave_sup_init() {
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    
    v_inputs_m_v_sleave_sup_init_wrap(&vlidort_sup_lcl);
    

  }

  void ss_sup_init() {
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    
    v_inputs_m_v_ss_sup_init_wrap(&vlidort_sup_lcl);
    

  }

  void sup_init() {
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    
    v_inputs_m_v_sup_init_wrap(&vlidort_sup_lcl);
    

  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_inputstatus_lcl = vlidort_inputstatus_->fortran_type_ptr();

    v_inputs_m_read_wrap(filename_lcl, &filename_in_len, &vlidort_sup_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_inputstatus_lcl);
    
  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) const {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_inputstatus_lcl = vlidort_inputstatus_->fortran_type_ptr();

    v_inputs_m_write_wrap(filename_lcl, &filename_in_len, &vlidort_sup_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_inputstatus_lcl);
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "VLidort_Inputs:" << std::endl
      << "        vlidort_sup: " << vlidort_sup()  << std::endl
      << "      vlidort_fixin: " << vlidort_fixin()  << std::endl
      << "      vlidort_modin: " << vlidort_modin()  << std::endl
      << "vlidort_inputstatus: " << vlidort_inputstatus()  << std::endl;

  }

private:
  boost::shared_ptr<VLidort_Sup_Inout> vlidort_sup_;
  boost::shared_ptr<VLidort_Fixed_Inputs> vlidort_fixin_;
  boost::shared_ptr<VLidort_Modified_Inputs> vlidort_modin_;
  boost::shared_ptr<VLidort_Input_Exception_Handling> vlidort_inputstatus_;

  // Serialization support
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "vlidort_masters_m" in file: "vlidort_masters.f90"
//-----------------------------------------------------------------------

extern "C" {
  void v_masters_m_read_wrap(const char* filename, const int* filename_len, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in);
  void v_masters_m_write_wrap(const char* filename, const int* filename_len, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in);
  void v_masters_m_v_master_wrap(const bool* do_debug_input_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in);
}

class VLidort_Masters : public virtual Spurr_Masters_Base {

public:
  VLidort_Masters() 
  { 
    
    // Initialize type pointers
    vlidort_fixin_.reset( new VLidort_Fixed_Inputs() );
    vlidort_modin_.reset( new VLidort_Modified_Inputs() );
    vlidort_sup_.reset( new VLidort_Sup_Inout() );
    vlidort_out_.reset( new VLidort_Outputs() );
    
  }

  virtual ~VLidort_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  VLidort_Fixed_Inputs& vlidort_fixin() {
    return *vlidort_fixin_;
  }

  const VLidort_Fixed_Inputs& vlidort_fixin() const {
    return *vlidort_fixin_;
  }

  Spurr_Fixed_Inputs_Base& fixed_inputs_base() {
    return *vlidort_fixin_;
  }

  const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const {
    return *vlidort_fixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Inputs>& vlidort_fixin_ptr() {
    return vlidort_fixin_;
  }

  void vlidort_fixin(VLidort_Fixed_Inputs& vlidort_fixin_in) {
    void* src_ptr = vlidort_fixin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_fixin_->fortran_type_ptr();
    vlidort_fixed_inputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Modified_Inputs& vlidort_modin() {
    return *vlidort_modin_;
  }

  const VLidort_Modified_Inputs& vlidort_modin() const {
    return *vlidort_modin_;
  }

  Spurr_Modified_Inputs_Base& modified_inputs_base() {
    return *vlidort_modin_;
  }

  const Spurr_Modified_Inputs_Base& modified_inputs_base() const {
    return *vlidort_modin_;
  }

  boost::shared_ptr<VLidort_Modified_Inputs>& vlidort_modin_ptr() {
    return vlidort_modin_;
  }

  void vlidort_modin(VLidort_Modified_Inputs& vlidort_modin_in) {
    void* src_ptr = vlidort_modin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_modin_->fortran_type_ptr();
    vlidort_modified_inputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Sup_Inout& vlidort_sup() {
    return *vlidort_sup_;
  }

  const VLidort_Sup_Inout& vlidort_sup() const {
    return *vlidort_sup_;
  }

  Spurr_Sup_Inout_Base& sup_inout_base() {
    return *vlidort_sup_;
  }

  const Spurr_Sup_Inout_Base& sup_inout_base() const {
    return *vlidort_sup_;
  }

  boost::shared_ptr<VLidort_Sup_Inout>& vlidort_sup_ptr() {
    return vlidort_sup_;
  }

  void vlidort_sup(VLidort_Sup_Inout& vlidort_sup_in) {
    void* src_ptr = vlidort_sup_in.fortran_type_ptr();
    void* dst_ptr = vlidort_sup_->fortran_type_ptr();
    vlidort_sup_inout_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Outputs& vlidort_out() {
    return *vlidort_out_;
  }

  const VLidort_Outputs& vlidort_out() const {
    return *vlidort_out_;
  }

  Spurr_Outputs_Base& outputs_base() {
    return *vlidort_out_;
  }

  const Spurr_Outputs_Base& outputs_base() const {
    return *vlidort_out_;
  }

  boost::shared_ptr<VLidort_Outputs>& vlidort_out_ptr() {
    return vlidort_out_;
  }

  
  void run(const bool& do_debug_input_in) {
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();
    
    v_masters_m_v_master_wrap(&do_debug_input_in, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl);
    
    VLidort_Pars vlid_pars = VLidort_Pars::instance();
    if( vlidort_out().status().ts_status_inputcheck() != vlid_pars.vlidort_success() ||
        vlidort_out().status().ts_status_calculation() != vlid_pars.vlidort_success() ) {
       std::stringstream err_msg;
       err_msg << "VLIDORT Error at " << __FILE__ << ":" << __LINE__ << std::endl;
       // Output the full details of the error message to stderr since the exception
       // class may truncate the message
       std::cerr << err_msg.str();
       std::cerr << vlidort_out().status();
       throw Exception(err_msg.str());
    }

  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();

    v_masters_m_read_wrap(filename_lcl, &filename_in_len, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl);
    
  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) const {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();

    v_masters_m_write_wrap(filename_lcl, &filename_in_len, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl);
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "VLidort_Masters:" << std::endl
      << "vlidort_fixin: " << vlidort_fixin()  << std::endl
      << "vlidort_modin: " << vlidort_modin()  << std::endl
      << "  vlidort_sup: " << vlidort_sup()  << std::endl
      << "  vlidort_out: " << vlidort_out()  << std::endl;

  }

private:
  boost::shared_ptr<VLidort_Fixed_Inputs> vlidort_fixin_;
  boost::shared_ptr<VLidort_Modified_Inputs> vlidort_modin_;
  boost::shared_ptr<VLidort_Sup_Inout> vlidort_sup_;
  boost::shared_ptr<VLidort_Outputs> vlidort_out_;

  // Serialization support
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "vlidort_l_inputs_m" in file: "vlidort_l_inputs.f90"
//-----------------------------------------------------------------------

extern "C" {
  void v_l_inputs_m_read_wrap(const char* filename, const int* filename_len, void** vlidort_linsup_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_inputstatus_in);
  void v_l_inputs_m_write_wrap(const char* filename, const int* filename_len, void** vlidort_linsup_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_inputstatus_in);
  void v_l_inputs_m_v_brdf_linsup_init_wrap(void** vlidort_linsup_in);
  void v_l_inputs_m_v_l_input_master_wrap(const int* filnam_in_len, const char* filnam_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_inputstatus_in);
  void v_l_inputs_m_v_linsup_init_wrap(void** vlidort_linsup_in);
  void v_l_inputs_m_v_sleave_linsup_init_wrap(void** vlidort_linsup_in);
  void v_l_inputs_m_v_ss_linsup_init_wrap(void** vlidort_linsup_in);
}

class VLidort_L_Inputs : public virtual Spurr_L_Inputs_Base {

public:
  VLidort_L_Inputs() 
  { 
    
    // Initialize type pointers
    vlidort_linsup_.reset( new VLidort_Linsup_Inout() );
    vlidort_fixin_.reset( new VLidort_Fixed_Inputs() );
    vlidort_modin_.reset( new VLidort_Modified_Inputs() );
    vlidort_linfixin_.reset( new VLidort_Fixed_Lininputs() );
    vlidort_linmodin_.reset( new VLidort_Modified_Lininputs() );
    vlidort_inputstatus_.reset( new VLidort_Input_Exception_Handling() );
    
  }

  virtual ~VLidort_L_Inputs() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  VLidort_Linsup_Inout& vlidort_linsup() {
    return *vlidort_linsup_;
  }

  const VLidort_Linsup_Inout& vlidort_linsup() const {
    return *vlidort_linsup_;
  }

  Spurr_Lin_Sup_Inout_Base& linsup_inout_base() {
    return *vlidort_linsup_;
  }

  const Spurr_Lin_Sup_Inout_Base& linsup_inout_base() const {
    return *vlidort_linsup_;
  }

  boost::shared_ptr<VLidort_Linsup_Inout>& vlidort_linsup_ptr() {
    return vlidort_linsup_;
  }

  void vlidort_linsup(VLidort_Linsup_Inout& vlidort_linsup_in) {
    void* src_ptr = vlidort_linsup_in.fortran_type_ptr();
    void* dst_ptr = vlidort_linsup_->fortran_type_ptr();
    vlidort_linsup_inout_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Fixed_Inputs& vlidort_fixin() {
    return *vlidort_fixin_;
  }

  const VLidort_Fixed_Inputs& vlidort_fixin() const {
    return *vlidort_fixin_;
  }

  Spurr_Fixed_Inputs_Base& fixed_inputs_base() {
    return *vlidort_fixin_;
  }

  const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const {
    return *vlidort_fixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Inputs>& vlidort_fixin_ptr() {
    return vlidort_fixin_;
  }

  
  VLidort_Modified_Inputs& vlidort_modin() {
    return *vlidort_modin_;
  }

  const VLidort_Modified_Inputs& vlidort_modin() const {
    return *vlidort_modin_;
  }

  Spurr_Modified_Inputs_Base& modified_inputs_base() {
    return *vlidort_modin_;
  }

  const Spurr_Modified_Inputs_Base& modified_inputs_base() const {
    return *vlidort_modin_;
  }

  boost::shared_ptr<VLidort_Modified_Inputs>& vlidort_modin_ptr() {
    return vlidort_modin_;
  }

  
  VLidort_Fixed_Lininputs& vlidort_linfixin() {
    return *vlidort_linfixin_;
  }

  const VLidort_Fixed_Lininputs& vlidort_linfixin() const {
    return *vlidort_linfixin_;
  }

  Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() {
    return *vlidort_linfixin_;
  }

  const Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() const {
    return *vlidort_linfixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Lininputs>& vlidort_linfixin_ptr() {
    return vlidort_linfixin_;
  }

  
  VLidort_Modified_Lininputs& vlidort_linmodin() {
    return *vlidort_linmodin_;
  }

  const VLidort_Modified_Lininputs& vlidort_linmodin() const {
    return *vlidort_linmodin_;
  }

  Spurr_Modified_Lininputs_Base& modified_lininputs_base() {
    return *vlidort_linmodin_;
  }

  const Spurr_Modified_Lininputs_Base& modified_lininputs_base() const {
    return *vlidort_linmodin_;
  }

  boost::shared_ptr<VLidort_Modified_Lininputs>& vlidort_linmodin_ptr() {
    return vlidort_linmodin_;
  }

  
  VLidort_Input_Exception_Handling& vlidort_inputstatus() {
    return *vlidort_inputstatus_;
  }

  const VLidort_Input_Exception_Handling& vlidort_inputstatus() const {
    return *vlidort_inputstatus_;
  }

  Spurr_Input_Exception_Handling_Base& input_exception_handling_base() {
    return *vlidort_inputstatus_;
  }

  const Spurr_Input_Exception_Handling_Base& input_exception_handling_base() const {
    return *vlidort_inputstatus_;
  }

  boost::shared_ptr<VLidort_Input_Exception_Handling>& vlidort_inputstatus_ptr() {
    return vlidort_inputstatus_;
  }

  
  void brdf_linsup_init() {
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    
    v_l_inputs_m_v_brdf_linsup_init_wrap(&vlidort_linsup_lcl);
    

  }

  void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_inputstatus_lcl = vlidort_inputstatus_->fortran_type_ptr();
    
    v_l_inputs_m_v_l_input_master_wrap(&filnam_in_len, filnam_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_inputstatus_lcl);
    

  }

  void linsup_init() {
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    
    v_l_inputs_m_v_linsup_init_wrap(&vlidort_linsup_lcl);
    

  }

  void sleave_linsup_init() {
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    
    v_l_inputs_m_v_sleave_linsup_init_wrap(&vlidort_linsup_lcl);
    

  }

  void ss_linsup_init() {
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    
    v_l_inputs_m_v_ss_linsup_init_wrap(&vlidort_linsup_lcl);
    

  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_inputstatus_lcl = vlidort_inputstatus_->fortran_type_ptr();

    v_l_inputs_m_read_wrap(filename_lcl, &filename_in_len, &vlidort_linsup_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_inputstatus_lcl);
    
  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) const {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_inputstatus_lcl = vlidort_inputstatus_->fortran_type_ptr();

    v_l_inputs_m_write_wrap(filename_lcl, &filename_in_len, &vlidort_linsup_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_inputstatus_lcl);
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "VLidort_L_Inputs:" << std::endl
      << "     vlidort_linsup: " << vlidort_linsup()  << std::endl
      << "      vlidort_fixin: " << vlidort_fixin()  << std::endl
      << "      vlidort_modin: " << vlidort_modin()  << std::endl
      << "   vlidort_linfixin: " << vlidort_linfixin()  << std::endl
      << "   vlidort_linmodin: " << vlidort_linmodin()  << std::endl
      << "vlidort_inputstatus: " << vlidort_inputstatus()  << std::endl;

  }

private:
  boost::shared_ptr<VLidort_Linsup_Inout> vlidort_linsup_;
  boost::shared_ptr<VLidort_Fixed_Inputs> vlidort_fixin_;
  boost::shared_ptr<VLidort_Modified_Inputs> vlidort_modin_;
  boost::shared_ptr<VLidort_Fixed_Lininputs> vlidort_linfixin_;
  boost::shared_ptr<VLidort_Modified_Lininputs> vlidort_linmodin_;
  boost::shared_ptr<VLidort_Input_Exception_Handling> vlidort_inputstatus_;

  // Serialization support
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "vlidort_lcs_masters_m" in file: "vlidort_lcs_masters.f90"
//-----------------------------------------------------------------------

extern "C" {
  void v_lcs_masters_m_read_wrap(const char* filename, const int* filename_len, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_linsup_in, void** vlidort_linout_in);
  void v_lcs_masters_m_write_wrap(const char* filename, const int* filename_len, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_linsup_in, void** vlidort_linout_in);
  void v_lcs_masters_m_v_lcs_master_wrap(const bool* do_debug_input_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_linsup_in, void** vlidort_linout_in);
}

class VLidort_Lcs_Masters : public virtual Spurr_Lcs_Masters_Base {

public:
  VLidort_Lcs_Masters() 
  { 
    
    // Initialize type pointers
    vlidort_fixin_.reset( new VLidort_Fixed_Inputs() );
    vlidort_modin_.reset( new VLidort_Modified_Inputs() );
    vlidort_sup_.reset( new VLidort_Sup_Inout() );
    vlidort_out_.reset( new VLidort_Outputs() );
    vlidort_linfixin_.reset( new VLidort_Fixed_Lininputs() );
    vlidort_linmodin_.reset( new VLidort_Modified_Lininputs() );
    vlidort_linsup_.reset( new VLidort_Linsup_Inout() );
    vlidort_linout_.reset( new VLidort_Linoutputs() );
    
  }

  virtual ~VLidort_Lcs_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  VLidort_Fixed_Inputs& vlidort_fixin() {
    return *vlidort_fixin_;
  }

  const VLidort_Fixed_Inputs& vlidort_fixin() const {
    return *vlidort_fixin_;
  }

  Spurr_Fixed_Inputs_Base& fixed_inputs_base() {
    return *vlidort_fixin_;
  }

  const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const {
    return *vlidort_fixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Inputs>& vlidort_fixin_ptr() {
    return vlidort_fixin_;
  }

  void vlidort_fixin(VLidort_Fixed_Inputs& vlidort_fixin_in) {
    void* src_ptr = vlidort_fixin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_fixin_->fortran_type_ptr();
    vlidort_fixed_inputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Modified_Inputs& vlidort_modin() {
    return *vlidort_modin_;
  }

  const VLidort_Modified_Inputs& vlidort_modin() const {
    return *vlidort_modin_;
  }

  Spurr_Modified_Inputs_Base& modified_inputs_base() {
    return *vlidort_modin_;
  }

  const Spurr_Modified_Inputs_Base& modified_inputs_base() const {
    return *vlidort_modin_;
  }

  boost::shared_ptr<VLidort_Modified_Inputs>& vlidort_modin_ptr() {
    return vlidort_modin_;
  }

  void vlidort_modin(VLidort_Modified_Inputs& vlidort_modin_in) {
    void* src_ptr = vlidort_modin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_modin_->fortran_type_ptr();
    vlidort_modified_inputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Sup_Inout& vlidort_sup() {
    return *vlidort_sup_;
  }

  const VLidort_Sup_Inout& vlidort_sup() const {
    return *vlidort_sup_;
  }

  Spurr_Sup_Inout_Base& sup_inout_base() {
    return *vlidort_sup_;
  }

  const Spurr_Sup_Inout_Base& sup_inout_base() const {
    return *vlidort_sup_;
  }

  boost::shared_ptr<VLidort_Sup_Inout>& vlidort_sup_ptr() {
    return vlidort_sup_;
  }

  void vlidort_sup(VLidort_Sup_Inout& vlidort_sup_in) {
    void* src_ptr = vlidort_sup_in.fortran_type_ptr();
    void* dst_ptr = vlidort_sup_->fortran_type_ptr();
    vlidort_sup_inout_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Outputs& vlidort_out() {
    return *vlidort_out_;
  }

  const VLidort_Outputs& vlidort_out() const {
    return *vlidort_out_;
  }

  Spurr_Outputs_Base& outputs_base() {
    return *vlidort_out_;
  }

  const Spurr_Outputs_Base& outputs_base() const {
    return *vlidort_out_;
  }

  boost::shared_ptr<VLidort_Outputs>& vlidort_out_ptr() {
    return vlidort_out_;
  }

  
  VLidort_Fixed_Lininputs& vlidort_linfixin() {
    return *vlidort_linfixin_;
  }

  const VLidort_Fixed_Lininputs& vlidort_linfixin() const {
    return *vlidort_linfixin_;
  }

  Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() {
    return *vlidort_linfixin_;
  }

  const Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() const {
    return *vlidort_linfixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Lininputs>& vlidort_linfixin_ptr() {
    return vlidort_linfixin_;
  }

  void vlidort_linfixin(VLidort_Fixed_Lininputs& vlidort_linfixin_in) {
    void* src_ptr = vlidort_linfixin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_linfixin_->fortran_type_ptr();
    vlidort_fixed_lininputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Modified_Lininputs& vlidort_linmodin() {
    return *vlidort_linmodin_;
  }

  const VLidort_Modified_Lininputs& vlidort_linmodin() const {
    return *vlidort_linmodin_;
  }

  Spurr_Modified_Lininputs_Base& modified_lininputs_base() {
    return *vlidort_linmodin_;
  }

  const Spurr_Modified_Lininputs_Base& modified_lininputs_base() const {
    return *vlidort_linmodin_;
  }

  boost::shared_ptr<VLidort_Modified_Lininputs>& vlidort_linmodin_ptr() {
    return vlidort_linmodin_;
  }

  void vlidort_linmodin(VLidort_Modified_Lininputs& vlidort_linmodin_in) {
    void* src_ptr = vlidort_linmodin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_linmodin_->fortran_type_ptr();
    vlidort_modified_lininputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Linsup_Inout& vlidort_linsup() {
    return *vlidort_linsup_;
  }

  const VLidort_Linsup_Inout& vlidort_linsup() const {
    return *vlidort_linsup_;
  }

  Spurr_Lin_Sup_Inout_Base& linsup_inout_base() {
    return *vlidort_linsup_;
  }

  const Spurr_Lin_Sup_Inout_Base& linsup_inout_base() const {
    return *vlidort_linsup_;
  }

  boost::shared_ptr<VLidort_Linsup_Inout>& vlidort_linsup_ptr() {
    return vlidort_linsup_;
  }

  void vlidort_linsup(VLidort_Linsup_Inout& vlidort_linsup_in) {
    void* src_ptr = vlidort_linsup_in.fortran_type_ptr();
    void* dst_ptr = vlidort_linsup_->fortran_type_ptr();
    vlidort_linsup_inout_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Linoutputs& vlidort_linout() {
    return *vlidort_linout_;
  }

  const VLidort_Linoutputs& vlidort_linout() const {
    return *vlidort_linout_;
  }

  Spurr_Linoutputs_Base& linoutputs_base() {
    return *vlidort_linout_;
  }

  const Spurr_Linoutputs_Base& linoutputs_base() const {
    return *vlidort_linout_;
  }

  boost::shared_ptr<VLidort_Linoutputs>& vlidort_linout_ptr() {
    return vlidort_linout_;
  }

  
  void run(const bool& do_debug_input_in) {
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    void* vlidort_linout_lcl = vlidort_linout_->fortran_type_ptr();
    
    v_lcs_masters_m_v_lcs_master_wrap(&do_debug_input_in, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_linsup_lcl, &vlidort_linout_lcl);
    
    VLidort_Pars vlid_pars = VLidort_Pars::instance();
    if( vlidort_out().status().ts_status_inputcheck() != vlid_pars.vlidort_success() ||
        vlidort_out().status().ts_status_calculation() != vlid_pars.vlidort_success() ) {
       std::stringstream err_msg;
       err_msg << "VLIDORT Error at " << __FILE__ << ":" << __LINE__ << std::endl;
       // Output the full details of the error message to stderr since the exception
       // class may truncate the message
       std::cerr << err_msg.str();
       std::cerr << vlidort_out().status();
       throw Exception(err_msg.str());
    }

  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    void* vlidort_linout_lcl = vlidort_linout_->fortran_type_ptr();

    v_lcs_masters_m_read_wrap(filename_lcl, &filename_in_len, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_linsup_lcl, &vlidort_linout_lcl);
    
  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) const {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    void* vlidort_linout_lcl = vlidort_linout_->fortran_type_ptr();

    v_lcs_masters_m_write_wrap(filename_lcl, &filename_in_len, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_linsup_lcl, &vlidort_linout_lcl);
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "VLidort_Lcs_Masters:" << std::endl
      << "   vlidort_fixin: " << vlidort_fixin()  << std::endl
      << "   vlidort_modin: " << vlidort_modin()  << std::endl
      << "     vlidort_sup: " << vlidort_sup()  << std::endl
      << "     vlidort_out: " << vlidort_out()  << std::endl
      << "vlidort_linfixin: " << vlidort_linfixin()  << std::endl
      << "vlidort_linmodin: " << vlidort_linmodin()  << std::endl
      << "  vlidort_linsup: " << vlidort_linsup()  << std::endl
      << "  vlidort_linout: " << vlidort_linout()  << std::endl;

  }

private:
  boost::shared_ptr<VLidort_Fixed_Inputs> vlidort_fixin_;
  boost::shared_ptr<VLidort_Modified_Inputs> vlidort_modin_;
  boost::shared_ptr<VLidort_Sup_Inout> vlidort_sup_;
  boost::shared_ptr<VLidort_Outputs> vlidort_out_;
  boost::shared_ptr<VLidort_Fixed_Lininputs> vlidort_linfixin_;
  boost::shared_ptr<VLidort_Modified_Lininputs> vlidort_linmodin_;
  boost::shared_ptr<VLidort_Linsup_Inout> vlidort_linsup_;
  boost::shared_ptr<VLidort_Linoutputs> vlidort_linout_;

  // Serialization support
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "vlidort_lps_masters_m" in file: "vlidort_lps_masters.f90"
//-----------------------------------------------------------------------

extern "C" {
  void v_lps_masters_m_read_wrap(const char* filename, const int* filename_len, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_linsup_in, void** vlidort_linout_in);
  void v_lps_masters_m_write_wrap(const char* filename, const int* filename_len, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_linsup_in, void** vlidort_linout_in);
  void v_lps_masters_m_v_lps_master_wrap(const bool* do_debug_input_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vlidort_out_in, void** vlidort_linfixin_in, void** vlidort_linmodin_in, void** vlidort_linsup_in, void** vlidort_linout_in);
}

class VLidort_Lps_Masters : public virtual Spurr_Lps_Masters_Base {

public:
  VLidort_Lps_Masters() 
  { 
    
    // Initialize type pointers
    vlidort_fixin_.reset( new VLidort_Fixed_Inputs() );
    vlidort_modin_.reset( new VLidort_Modified_Inputs() );
    vlidort_sup_.reset( new VLidort_Sup_Inout() );
    vlidort_out_.reset( new VLidort_Outputs() );
    vlidort_linfixin_.reset( new VLidort_Fixed_Lininputs() );
    vlidort_linmodin_.reset( new VLidort_Modified_Lininputs() );
    vlidort_linsup_.reset( new VLidort_Linsup_Inout() );
    vlidort_linout_.reset( new VLidort_Linoutputs() );
    
  }

  virtual ~VLidort_Lps_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  VLidort_Fixed_Inputs& vlidort_fixin() {
    return *vlidort_fixin_;
  }

  const VLidort_Fixed_Inputs& vlidort_fixin() const {
    return *vlidort_fixin_;
  }

  Spurr_Fixed_Inputs_Base& fixed_inputs_base() {
    return *vlidort_fixin_;
  }

  const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const {
    return *vlidort_fixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Inputs>& vlidort_fixin_ptr() {
    return vlidort_fixin_;
  }

  void vlidort_fixin(VLidort_Fixed_Inputs& vlidort_fixin_in) {
    void* src_ptr = vlidort_fixin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_fixin_->fortran_type_ptr();
    vlidort_fixed_inputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Modified_Inputs& vlidort_modin() {
    return *vlidort_modin_;
  }

  const VLidort_Modified_Inputs& vlidort_modin() const {
    return *vlidort_modin_;
  }

  Spurr_Modified_Inputs_Base& modified_inputs_base() {
    return *vlidort_modin_;
  }

  const Spurr_Modified_Inputs_Base& modified_inputs_base() const {
    return *vlidort_modin_;
  }

  boost::shared_ptr<VLidort_Modified_Inputs>& vlidort_modin_ptr() {
    return vlidort_modin_;
  }

  void vlidort_modin(VLidort_Modified_Inputs& vlidort_modin_in) {
    void* src_ptr = vlidort_modin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_modin_->fortran_type_ptr();
    vlidort_modified_inputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Sup_Inout& vlidort_sup() {
    return *vlidort_sup_;
  }

  const VLidort_Sup_Inout& vlidort_sup() const {
    return *vlidort_sup_;
  }

  Spurr_Sup_Inout_Base& sup_inout_base() {
    return *vlidort_sup_;
  }

  const Spurr_Sup_Inout_Base& sup_inout_base() const {
    return *vlidort_sup_;
  }

  boost::shared_ptr<VLidort_Sup_Inout>& vlidort_sup_ptr() {
    return vlidort_sup_;
  }

  void vlidort_sup(VLidort_Sup_Inout& vlidort_sup_in) {
    void* src_ptr = vlidort_sup_in.fortran_type_ptr();
    void* dst_ptr = vlidort_sup_->fortran_type_ptr();
    vlidort_sup_inout_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Outputs& vlidort_out() {
    return *vlidort_out_;
  }

  const VLidort_Outputs& vlidort_out() const {
    return *vlidort_out_;
  }

  Spurr_Outputs_Base& outputs_base() {
    return *vlidort_out_;
  }

  const Spurr_Outputs_Base& outputs_base() const {
    return *vlidort_out_;
  }

  boost::shared_ptr<VLidort_Outputs>& vlidort_out_ptr() {
    return vlidort_out_;
  }

  
  VLidort_Fixed_Lininputs& vlidort_linfixin() {
    return *vlidort_linfixin_;
  }

  const VLidort_Fixed_Lininputs& vlidort_linfixin() const {
    return *vlidort_linfixin_;
  }

  Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() {
    return *vlidort_linfixin_;
  }

  const Spurr_Fixed_Lininputs_Base& fixed_lininputs_base() const {
    return *vlidort_linfixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Lininputs>& vlidort_linfixin_ptr() {
    return vlidort_linfixin_;
  }

  void vlidort_linfixin(VLidort_Fixed_Lininputs& vlidort_linfixin_in) {
    void* src_ptr = vlidort_linfixin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_linfixin_->fortran_type_ptr();
    vlidort_fixed_lininputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Modified_Lininputs& vlidort_linmodin() {
    return *vlidort_linmodin_;
  }

  const VLidort_Modified_Lininputs& vlidort_linmodin() const {
    return *vlidort_linmodin_;
  }

  Spurr_Modified_Lininputs_Base& modified_lininputs_base() {
    return *vlidort_linmodin_;
  }

  const Spurr_Modified_Lininputs_Base& modified_lininputs_base() const {
    return *vlidort_linmodin_;
  }

  boost::shared_ptr<VLidort_Modified_Lininputs>& vlidort_linmodin_ptr() {
    return vlidort_linmodin_;
  }

  void vlidort_linmodin(VLidort_Modified_Lininputs& vlidort_linmodin_in) {
    void* src_ptr = vlidort_linmodin_in.fortran_type_ptr();
    void* dst_ptr = vlidort_linmodin_->fortran_type_ptr();
    vlidort_modified_lininputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Linsup_Inout& vlidort_linsup() {
    return *vlidort_linsup_;
  }

  const VLidort_Linsup_Inout& vlidort_linsup() const {
    return *vlidort_linsup_;
  }

  Spurr_Lin_Sup_Inout_Base& linsup_inout_base() {
    return *vlidort_linsup_;
  }

  const Spurr_Lin_Sup_Inout_Base& linsup_inout_base() const {
    return *vlidort_linsup_;
  }

  boost::shared_ptr<VLidort_Linsup_Inout>& vlidort_linsup_ptr() {
    return vlidort_linsup_;
  }

  void vlidort_linsup(VLidort_Linsup_Inout& vlidort_linsup_in) {
    void* src_ptr = vlidort_linsup_in.fortran_type_ptr();
    void* dst_ptr = vlidort_linsup_->fortran_type_ptr();
    vlidort_linsup_inout_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Linoutputs& vlidort_linout() {
    return *vlidort_linout_;
  }

  const VLidort_Linoutputs& vlidort_linout() const {
    return *vlidort_linout_;
  }

  Spurr_Linoutputs_Base& linoutputs_base() {
    return *vlidort_linout_;
  }

  const Spurr_Linoutputs_Base& linoutputs_base() const {
    return *vlidort_linout_;
  }

  boost::shared_ptr<VLidort_Linoutputs>& vlidort_linout_ptr() {
    return vlidort_linout_;
  }

  
  void run(const bool& do_debug_input_in) {
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    void* vlidort_linout_lcl = vlidort_linout_->fortran_type_ptr();
    
    v_lps_masters_m_v_lps_master_wrap(&do_debug_input_in, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_linsup_lcl, &vlidort_linout_lcl);
    
    VLidort_Pars vlid_pars = VLidort_Pars::instance();
    if( vlidort_out().status().ts_status_inputcheck() != vlid_pars.vlidort_success() ||
        vlidort_out().status().ts_status_calculation() != vlid_pars.vlidort_success() ) {
       std::stringstream err_msg;
       err_msg << "VLIDORT Error at " << __FILE__ << ":" << __LINE__ << std::endl;
       // Output the full details of the error message to stderr since the exception
       // class may truncate the message
       std::cerr << err_msg.str();
       std::cerr << vlidort_out().status();
       throw Exception(err_msg.str());
    }

  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    void* vlidort_linout_lcl = vlidort_linout_->fortran_type_ptr();

    v_lps_masters_m_read_wrap(filename_lcl, &filename_in_len, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_linsup_lcl, &vlidort_linout_lcl);
    
  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) const {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vlidort_out_lcl = vlidort_out_->fortran_type_ptr();
    void* vlidort_linfixin_lcl = vlidort_linfixin_->fortran_type_ptr();
    void* vlidort_linmodin_lcl = vlidort_linmodin_->fortran_type_ptr();
    void* vlidort_linsup_lcl = vlidort_linsup_->fortran_type_ptr();
    void* vlidort_linout_lcl = vlidort_linout_->fortran_type_ptr();

    v_lps_masters_m_write_wrap(filename_lcl, &filename_in_len, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vlidort_out_lcl, &vlidort_linfixin_lcl, &vlidort_linmodin_lcl, &vlidort_linsup_lcl, &vlidort_linout_lcl);
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "VLidort_Lps_Masters:" << std::endl
      << "   vlidort_fixin: " << vlidort_fixin()  << std::endl
      << "   vlidort_modin: " << vlidort_modin()  << std::endl
      << "     vlidort_sup: " << vlidort_sup()  << std::endl
      << "     vlidort_out: " << vlidort_out()  << std::endl
      << "vlidort_linfixin: " << vlidort_linfixin()  << std::endl
      << "vlidort_linmodin: " << vlidort_linmodin()  << std::endl
      << "  vlidort_linsup: " << vlidort_linsup()  << std::endl
      << "  vlidort_linout: " << vlidort_linout()  << std::endl;

  }

private:
  boost::shared_ptr<VLidort_Fixed_Inputs> vlidort_fixin_;
  boost::shared_ptr<VLidort_Modified_Inputs> vlidort_modin_;
  boost::shared_ptr<VLidort_Sup_Inout> vlidort_sup_;
  boost::shared_ptr<VLidort_Outputs> vlidort_out_;
  boost::shared_ptr<VLidort_Fixed_Lininputs> vlidort_linfixin_;
  boost::shared_ptr<VLidort_Modified_Lininputs> vlidort_linmodin_;
  boost::shared_ptr<VLidort_Linsup_Inout> vlidort_linsup_;
  boost::shared_ptr<VLidort_Linoutputs> vlidort_linout_;

  // Serialization support
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "vlidort_vbrdf_sup_accessories_m" in file: "vlidort_vbrdf_sup_accessories.f90"
//-----------------------------------------------------------------------

extern "C" {
  void v_vbrdf_sup_accessories_m_read_wrap(const char* filename, const int* filename_len, void** vbrdf_sup_out_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vbrdf_sup_in_in, void** vlidort_vbrdfcheck_status_in);
  void v_vbrdf_sup_accessories_m_write_wrap(const char* filename, const int* filename_len, void** vbrdf_sup_out_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in, void** vbrdf_sup_in_in, void** vlidort_vbrdfcheck_status_in);
  void v_vbrdf_sup_accessories_m_set_v_vbrdf_inputs_wrap(void** vbrdf_sup_out_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_sup_in);
  void v_vbrdf_sup_accessories_m_v_vbrdf_input_check_wrap(void** vbrdf_sup_in_in, void** vlidort_fixin_in, void** vlidort_modin_in, void** vlidort_vbrdfcheck_status_in);
  void v_vbrdf_sup_accessories_m_v_vbrdf_input_check_error_wrap(const int* errorfile_in_len, const char* errorfile_in, void** vlidort_vbrdfcheck_status_in);
}

class VLidort_Vbrdf_Sup_Accessories : public virtual Spurr_Brdf_Sup_Accessories_Base {

public:
  VLidort_Vbrdf_Sup_Accessories(boost::shared_ptr<VLidort_Fixed_Inputs>& vlidort_fixin_in, boost::shared_ptr<VLidort_Modified_Inputs>& vlidort_modin_in, boost::shared_ptr<VBrdf_Sup_Inputs>& vbrdf_sup_in_in) : vlidort_fixin_(vlidort_fixin_in), vlidort_modin_(vlidort_modin_in), vbrdf_sup_in_(vbrdf_sup_in_in) 
  { 
    
    // Initialize type pointers
    vbrdf_sup_out_.reset( new VBrdf_Sup_Outputs() );
    vlidort_sup_.reset( new VLidort_Sup_Inout() );
    vlidort_vbrdfcheck_status_.reset( new VLidort_Exception_Handling() );
    
  }

  virtual ~VLidort_Vbrdf_Sup_Accessories() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  VBrdf_Sup_Outputs& vbrdf_sup_out() {
    return *vbrdf_sup_out_;
  }

  const VBrdf_Sup_Outputs& vbrdf_sup_out() const {
    return *vbrdf_sup_out_;
  }

  Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() {
    return *vbrdf_sup_out_;
  }

  const Spurr_Brdf_Sup_Outputs_Base& brdf_sup_outputs_base() const {
    return *vbrdf_sup_out_;
  }

  boost::shared_ptr<VBrdf_Sup_Outputs>& vbrdf_sup_out_ptr() {
    return vbrdf_sup_out_;
  }

  void vbrdf_sup_out(VBrdf_Sup_Outputs& vbrdf_sup_out_in) {
    void* src_ptr = vbrdf_sup_out_in.fortran_type_ptr();
    void* dst_ptr = vbrdf_sup_out_->fortran_type_ptr();
    vbrdf_sup_outputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  VLidort_Fixed_Inputs& vlidort_fixin() {
    return *vlidort_fixin_;
  }

  const VLidort_Fixed_Inputs& vlidort_fixin() const {
    return *vlidort_fixin_;
  }

  Spurr_Fixed_Inputs_Base& fixed_inputs_base() {
    return *vlidort_fixin_;
  }

  const Spurr_Fixed_Inputs_Base& fixed_inputs_base() const {
    return *vlidort_fixin_;
  }

  boost::shared_ptr<VLidort_Fixed_Inputs>& vlidort_fixin_ptr() {
    return vlidort_fixin_;
  }

  
  VLidort_Modified_Inputs& vlidort_modin() {
    return *vlidort_modin_;
  }

  const VLidort_Modified_Inputs& vlidort_modin() const {
    return *vlidort_modin_;
  }

  Spurr_Modified_Inputs_Base& modified_inputs_base() {
    return *vlidort_modin_;
  }

  const Spurr_Modified_Inputs_Base& modified_inputs_base() const {
    return *vlidort_modin_;
  }

  boost::shared_ptr<VLidort_Modified_Inputs>& vlidort_modin_ptr() {
    return vlidort_modin_;
  }

  
  VLidort_Sup_Inout& vlidort_sup() {
    return *vlidort_sup_;
  }

  const VLidort_Sup_Inout& vlidort_sup() const {
    return *vlidort_sup_;
  }

  Spurr_Sup_Inout_Base& sup_inout_base() {
    return *vlidort_sup_;
  }

  const Spurr_Sup_Inout_Base& sup_inout_base() const {
    return *vlidort_sup_;
  }

  boost::shared_ptr<VLidort_Sup_Inout>& vlidort_sup_ptr() {
    return vlidort_sup_;
  }

  void vlidort_sup(VLidort_Sup_Inout& vlidort_sup_in) {
    void* src_ptr = vlidort_sup_in.fortran_type_ptr();
    void* dst_ptr = vlidort_sup_->fortran_type_ptr();
    vlidort_sup_inout_c_copy(&src_ptr, &dst_ptr);
  }

  
  VBrdf_Sup_Inputs& vbrdf_sup_in() {
    return *vbrdf_sup_in_;
  }

  const VBrdf_Sup_Inputs& vbrdf_sup_in() const {
    return *vbrdf_sup_in_;
  }

  Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() {
    return *vbrdf_sup_in_;
  }

  const Spurr_Brdf_Sup_Inputs_Base& brdf_sup_inputs_base() const {
    return *vbrdf_sup_in_;
  }

  boost::shared_ptr<VBrdf_Sup_Inputs>& vbrdf_sup_in_ptr() {
    return vbrdf_sup_in_;
  }

  
  VLidort_Exception_Handling& vlidort_vbrdfcheck_status() {
    return *vlidort_vbrdfcheck_status_;
  }

  const VLidort_Exception_Handling& vlidort_vbrdfcheck_status() const {
    return *vlidort_vbrdfcheck_status_;
  }

  Spurr_Exception_Handling_Base& exception_handling_base() {
    return *vlidort_vbrdfcheck_status_;
  }

  const Spurr_Exception_Handling_Base& exception_handling_base() const {
    return *vlidort_vbrdfcheck_status_;
  }

  boost::shared_ptr<VLidort_Exception_Handling>& vlidort_vbrdfcheck_status_ptr() {
    return vlidort_vbrdfcheck_status_;
  }

  
  void set_vbrdf_inputs() {
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    
    v_vbrdf_sup_accessories_m_set_v_vbrdf_inputs_wrap(&vbrdf_sup_out_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl);
    

  }

  void vbrdf_input_check() {
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_vbrdfcheck_status_lcl = vlidort_vbrdfcheck_status_->fortran_type_ptr();
    
    v_vbrdf_sup_accessories_m_v_vbrdf_input_check_wrap(&vbrdf_sup_in_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_vbrdfcheck_status_lcl);
    

  }

  void vbrdf_input_check_error(const std::string& errorfile_in) {
    const char* errorfile_lcl = errorfile_in.c_str();
    int errorfile_in_len = (int) errorfile_in.size();
    void* vlidort_vbrdfcheck_status_lcl = vlidort_vbrdfcheck_status_->fortran_type_ptr();
    
    v_vbrdf_sup_accessories_m_v_vbrdf_input_check_error_wrap(&errorfile_in_len, errorfile_lcl, &vlidort_vbrdfcheck_status_lcl);
    

  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vlidort_vbrdfcheck_status_lcl = vlidort_vbrdfcheck_status_->fortran_type_ptr();

    v_vbrdf_sup_accessories_m_read_wrap(filename_lcl, &filename_in_len, &vbrdf_sup_out_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vbrdf_sup_in_lcl, &vlidort_vbrdfcheck_status_lcl);
    
  }

  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) const {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* vbrdf_sup_out_lcl = vbrdf_sup_out_->fortran_type_ptr();
    void* vlidort_fixin_lcl = vlidort_fixin_->fortran_type_ptr();
    void* vlidort_modin_lcl = vlidort_modin_->fortran_type_ptr();
    void* vlidort_sup_lcl = vlidort_sup_->fortran_type_ptr();
    void* vbrdf_sup_in_lcl = vbrdf_sup_in_->fortran_type_ptr();
    void* vlidort_vbrdfcheck_status_lcl = vlidort_vbrdfcheck_status_->fortran_type_ptr();

    v_vbrdf_sup_accessories_m_write_wrap(filename_lcl, &filename_in_len, &vbrdf_sup_out_lcl, &vlidort_fixin_lcl, &vlidort_modin_lcl, &vlidort_sup_lcl, &vbrdf_sup_in_lcl, &vlidort_vbrdfcheck_status_lcl);
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "VLidort_Vbrdf_Sup_Accessories:" << std::endl
      << "            vbrdf_sup_out: " << vbrdf_sup_out()  << std::endl
      << "            vlidort_fixin: " << vlidort_fixin()  << std::endl
      << "            vlidort_modin: " << vlidort_modin()  << std::endl
      << "              vlidort_sup: " << vlidort_sup()  << std::endl
      << "             vbrdf_sup_in: " << vbrdf_sup_in()  << std::endl
      << "vlidort_vbrdfcheck_status: " << vlidort_vbrdfcheck_status()  << std::endl;

  }

private:
  boost::shared_ptr<VBrdf_Sup_Outputs> vbrdf_sup_out_;
  boost::shared_ptr<VLidort_Fixed_Inputs> vlidort_fixin_;
  boost::shared_ptr<VLidort_Modified_Inputs> vlidort_modin_;
  boost::shared_ptr<VLidort_Sup_Inout> vlidort_sup_;
  boost::shared_ptr<VBrdf_Sup_Inputs> vbrdf_sup_in_;
  boost::shared_ptr<VLidort_Exception_Handling> vlidort_vbrdfcheck_status_;

  // Serialization support
  VLidort_Vbrdf_Sup_Accessories() {}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

}

// Serialization support
FP_EXPORT_KEY(VBrdf_Linsup_Masters)
FP_EXPORT_KEY(VBrdf_Sup_Masters)
FP_EXPORT_KEY(VLidort_Inputs)
FP_EXPORT_KEY(VLidort_Masters)
FP_EXPORT_KEY(VLidort_L_Inputs)
FP_EXPORT_KEY(VLidort_Lcs_Masters)
FP_EXPORT_KEY(VLidort_Lps_Masters)
FP_EXPORT_KEY(VLidort_Vbrdf_Sup_Accessories)

#endif