#ifndef LIDORT_INTERFACE_TYPES_H
#define LIDORT_INTERFACE_TYPES_H

#include <iostream>
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "fp_exception.h"
#include "spurr_interface_types.h"


/* This file was auto-generated */

namespace FullPhysics {

extern "C" {
  void set_lidort_pars(struct lidort_pars_m *lidort_pars_struct_c);
}

/* This struct cannot inherit any properties from another struct/class which is a requirement of Fortran interoperability */

struct lidort_pars_m {

  double bigexp;
  int bpdfndvi_idx;
  int bpdfsoil_idx;
  int bpdfvegn_idx;
  int coxmunk_idx;
  double deg_to_rad;
  int dnidx;
  double eps3;
  double eps4;
  double eps5;
  double four;
  double half;
  int hapke_idx;
  double hopital_tolerance;
  int lambertian_idx;
  int lidense_idx;
  int lidort_dbgunit;
  int lidort_debug;
  int lidort_errunit;
  int lidort_funit;
  int lidort_info;
  int lidort_inunit;
  int lidort_resunit;
  int lidort_scenunit;
  int lidort_serious;
  int lidort_success;
  char lidort_version_number[5];
  int lidort_warning;
  int lisparse_idx;
  int max_allstrms;
  int max_allstrms_p1;
  int max_atmoswfs;
  int max_brdf_kernels;
  int max_brdf_parameters;
  int max_directions;
  int max_geometries;
  int max_messages;
  int max_msrs_muquad;
  int max_msrs_phiquad;
  int max_partlayers;
  int max_sleavewfs;
  int max_surfacewfs;
  double max_tau_qpath;
  double max_tau_spath;
  double max_tau_upath;
  int max_taylor_terms;
  int max_thermal_coeffs;
  int max_user_levels;
  int max_user_obsgeoms;
  int max_user_relazms;
  int max_user_streams;
  int maxbandtotal;
  int maxbeams;
  int maxbrdf_idx;
  int maxfinelayers;
  int maxfourier;
  int maxlayers;
  int maxmoments;
  int maxmoments_input;
  int maxsthalf_brdf;
  int maxstreams;
  int maxstreams_2;
  int maxstreams_brdf;
  int maxstreams_p1;
  int maxstreams_scaling;
  int maxtotal;
  double minus_one;
  double minus_two;
  int modfresnel_idx;
  int newcmglint_idx;
  double omega_smallnum;
  double one;
  double onep5;
  double pi2;
  double pi4;
  double pie;
  double pio2;
  double pio4;
  double quarter;
  int rahman_idx;
  int rossthick_idx;
  int rossthin_idx;
  int roujean_idx;
  int rtkhotspot_idx;
  double smallnum;
  int snowbrdf_idx;
  double taylor_small;
  double three;
  double two;
  int upidx;
  double zero;
  
};

class Lidort_Pars : public Spurr_Pars_Base {
public:

  const double& bigexp() const {
    return transfer_struct_c.bigexp;
  }

  const int& bpdfndvi_idx() const {
    return transfer_struct_c.bpdfndvi_idx;
  }

  const int& bpdfsoil_idx() const {
    return transfer_struct_c.bpdfsoil_idx;
  }

  const int& bpdfvegn_idx() const {
    return transfer_struct_c.bpdfvegn_idx;
  }

  const int& coxmunk_idx() const {
    return transfer_struct_c.coxmunk_idx;
  }

  const double& deg_to_rad() const {
    return transfer_struct_c.deg_to_rad;
  }

  const int& dnidx() const {
    return transfer_struct_c.dnidx;
  }

  const double& eps3() const {
    return transfer_struct_c.eps3;
  }

  const double& eps4() const {
    return transfer_struct_c.eps4;
  }

  const double& eps5() const {
    return transfer_struct_c.eps5;
  }

  const double& four() const {
    return transfer_struct_c.four;
  }

  const double& half() const {
    return transfer_struct_c.half;
  }

  const int& hapke_idx() const {
    return transfer_struct_c.hapke_idx;
  }

  const double& hopital_tolerance() const {
    return transfer_struct_c.hopital_tolerance;
  }

  const int& lambertian_idx() const {
    return transfer_struct_c.lambertian_idx;
  }

  const int& lidense_idx() const {
    return transfer_struct_c.lidense_idx;
  }

  const int& lidort_dbgunit() const {
    return transfer_struct_c.lidort_dbgunit;
  }

  const int& lidort_debug() const {
    return transfer_struct_c.lidort_debug;
  }

  const int& lidort_errunit() const {
    return transfer_struct_c.lidort_errunit;
  }

  const int& lidort_funit() const {
    return transfer_struct_c.lidort_funit;
  }

  const int& lidort_info() const {
    return transfer_struct_c.lidort_info;
  }

  const int& lidort_inunit() const {
    return transfer_struct_c.lidort_inunit;
  }

  const int& lidort_resunit() const {
    return transfer_struct_c.lidort_resunit;
  }

  const int& lidort_scenunit() const {
    return transfer_struct_c.lidort_scenunit;
  }

  const int& lidort_serious() const {
    return transfer_struct_c.lidort_serious;
  }

  const int& lidort_success() const {
    return transfer_struct_c.lidort_success;
  }

  const std::string lidort_version_number() const {
    return std::string(transfer_struct_c.lidort_version_number);
  }

  const int& lidort_warning() const {
    return transfer_struct_c.lidort_warning;
  }

  const int& lisparse_idx() const {
    return transfer_struct_c.lisparse_idx;
  }

  const int& max_allstrms() const {
    return transfer_struct_c.max_allstrms;
  }

  const int& max_allstrms_p1() const {
    return transfer_struct_c.max_allstrms_p1;
  }

  const int& max_atmoswfs() const {
    return transfer_struct_c.max_atmoswfs;
  }

  const int& max_brdf_kernels() const {
    return transfer_struct_c.max_brdf_kernels;
  }

  const int& max_brdf_parameters() const {
    return transfer_struct_c.max_brdf_parameters;
  }

  const int& max_directions() const {
    return transfer_struct_c.max_directions;
  }

  const int& max_geometries() const {
    return transfer_struct_c.max_geometries;
  }

  const int& max_messages() const {
    return transfer_struct_c.max_messages;
  }

  const int& max_msrs_muquad() const {
    return transfer_struct_c.max_msrs_muquad;
  }

  const int& max_msrs_phiquad() const {
    return transfer_struct_c.max_msrs_phiquad;
  }

  const int& max_partlayers() const {
    return transfer_struct_c.max_partlayers;
  }

  const int& max_sleavewfs() const {
    return transfer_struct_c.max_sleavewfs;
  }

  const int& max_surfacewfs() const {
    return transfer_struct_c.max_surfacewfs;
  }

  const double& max_tau_qpath() const {
    return transfer_struct_c.max_tau_qpath;
  }

  const double& max_tau_spath() const {
    return transfer_struct_c.max_tau_spath;
  }

  const double& max_tau_upath() const {
    return transfer_struct_c.max_tau_upath;
  }

  const int& max_taylor_terms() const {
    return transfer_struct_c.max_taylor_terms;
  }

  const int& max_thermal_coeffs() const {
    return transfer_struct_c.max_thermal_coeffs;
  }

  const int& max_user_levels() const {
    return transfer_struct_c.max_user_levels;
  }

  const int& max_user_obsgeoms() const {
    return transfer_struct_c.max_user_obsgeoms;
  }

  const int& max_user_relazms() const {
    return transfer_struct_c.max_user_relazms;
  }

  const int& max_user_streams() const {
    return transfer_struct_c.max_user_streams;
  }

  const int& maxbandtotal() const {
    return transfer_struct_c.maxbandtotal;
  }

  const int& maxbeams() const {
    return transfer_struct_c.maxbeams;
  }

  const int& maxbrdf_idx() const {
    return transfer_struct_c.maxbrdf_idx;
  }

  const int& maxfinelayers() const {
    return transfer_struct_c.maxfinelayers;
  }

  const int& maxfourier() const {
    return transfer_struct_c.maxfourier;
  }

  const int& maxlayers() const {
    return transfer_struct_c.maxlayers;
  }

  const int& maxmoments() const {
    return transfer_struct_c.maxmoments;
  }

  const int& maxmoments_input() const {
    return transfer_struct_c.maxmoments_input;
  }

  const int& maxsthalf_brdf() const {
    return transfer_struct_c.maxsthalf_brdf;
  }

  const int& maxstreams() const {
    return transfer_struct_c.maxstreams;
  }

  const int& maxstreams_2() const {
    return transfer_struct_c.maxstreams_2;
  }

  const int& maxstreams_brdf() const {
    return transfer_struct_c.maxstreams_brdf;
  }

  const int& maxstreams_p1() const {
    return transfer_struct_c.maxstreams_p1;
  }

  const int& maxstreams_scaling() const {
    return transfer_struct_c.maxstreams_scaling;
  }

  const int& maxtotal() const {
    return transfer_struct_c.maxtotal;
  }

  const double& minus_one() const {
    return transfer_struct_c.minus_one;
  }

  const double& minus_two() const {
    return transfer_struct_c.minus_two;
  }

  const int& modfresnel_idx() const {
    return transfer_struct_c.modfresnel_idx;
  }

  const int& newcmglint_idx() const {
    return transfer_struct_c.newcmglint_idx;
  }

  const double& omega_smallnum() const {
    return transfer_struct_c.omega_smallnum;
  }

  const double& one() const {
    return transfer_struct_c.one;
  }

  const double& onep5() const {
    return transfer_struct_c.onep5;
  }

  const double& pi2() const {
    return transfer_struct_c.pi2;
  }

  const double& pi4() const {
    return transfer_struct_c.pi4;
  }

  const double& pie() const {
    return transfer_struct_c.pie;
  }

  const double& pio2() const {
    return transfer_struct_c.pio2;
  }

  const double& pio4() const {
    return transfer_struct_c.pio4;
  }

  const double& quarter() const {
    return transfer_struct_c.quarter;
  }

  const int& rahman_idx() const {
    return transfer_struct_c.rahman_idx;
  }

  const int& rossthick_idx() const {
    return transfer_struct_c.rossthick_idx;
  }

  const int& rossthin_idx() const {
    return transfer_struct_c.rossthin_idx;
  }

  const int& roujean_idx() const {
    return transfer_struct_c.roujean_idx;
  }

  const int& rtkhotspot_idx() const {
    return transfer_struct_c.rtkhotspot_idx;
  }

  const double& smallnum() const {
    return transfer_struct_c.smallnum;
  }

  const int& snowbrdf_idx() const {
    return transfer_struct_c.snowbrdf_idx;
  }

  const double& taylor_small() const {
    return transfer_struct_c.taylor_small;
  }

  const double& three() const {
    return transfer_struct_c.three;
  }

  const double& two() const {
    return transfer_struct_c.two;
  }

  const int& upidx() const {
    return transfer_struct_c.upidx;
  }

  const double& zero() const {
    return transfer_struct_c.zero;
  }

  static Lidort_Pars& instance() {
    static Lidort_Pars obj;
    return obj;
  }
  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Pars:" << std::endl
      << "               bigexp: " << bigexp()  << std::endl
      << "         bpdfndvi_idx: " << bpdfndvi_idx()  << std::endl
      << "         bpdfsoil_idx: " << bpdfsoil_idx()  << std::endl
      << "         bpdfvegn_idx: " << bpdfvegn_idx()  << std::endl
      << "          coxmunk_idx: " << coxmunk_idx()  << std::endl
      << "           deg_to_rad: " << deg_to_rad()  << std::endl
      << "                dnidx: " << dnidx()  << std::endl
      << "                 eps3: " << eps3()  << std::endl
      << "                 eps4: " << eps4()  << std::endl
      << "                 eps5: " << eps5()  << std::endl
      << "                 four: " << four()  << std::endl
      << "                 half: " << half()  << std::endl
      << "            hapke_idx: " << hapke_idx()  << std::endl
      << "    hopital_tolerance: " << hopital_tolerance()  << std::endl
      << "       lambertian_idx: " << lambertian_idx()  << std::endl
      << "          lidense_idx: " << lidense_idx()  << std::endl
      << "       lidort_dbgunit: " << lidort_dbgunit()  << std::endl
      << "         lidort_debug: " << lidort_debug()  << std::endl
      << "       lidort_errunit: " << lidort_errunit()  << std::endl
      << "         lidort_funit: " << lidort_funit()  << std::endl
      << "          lidort_info: " << lidort_info()  << std::endl
      << "        lidort_inunit: " << lidort_inunit()  << std::endl
      << "       lidort_resunit: " << lidort_resunit()  << std::endl
      << "      lidort_scenunit: " << lidort_scenunit()  << std::endl
      << "       lidort_serious: " << lidort_serious()  << std::endl
      << "       lidort_success: " << lidort_success()  << std::endl
      << "lidort_version_number: " << "\"" << lidort_version_number() << "\"" << std::endl
      << "       lidort_warning: " << lidort_warning()  << std::endl
      << "         lisparse_idx: " << lisparse_idx()  << std::endl
      << "         max_allstrms: " << max_allstrms()  << std::endl
      << "      max_allstrms_p1: " << max_allstrms_p1()  << std::endl
      << "         max_atmoswfs: " << max_atmoswfs()  << std::endl
      << "     max_brdf_kernels: " << max_brdf_kernels()  << std::endl
      << "  max_brdf_parameters: " << max_brdf_parameters()  << std::endl
      << "       max_directions: " << max_directions()  << std::endl
      << "       max_geometries: " << max_geometries()  << std::endl
      << "         max_messages: " << max_messages()  << std::endl
      << "      max_msrs_muquad: " << max_msrs_muquad()  << std::endl
      << "     max_msrs_phiquad: " << max_msrs_phiquad()  << std::endl
      << "       max_partlayers: " << max_partlayers()  << std::endl
      << "        max_sleavewfs: " << max_sleavewfs()  << std::endl
      << "       max_surfacewfs: " << max_surfacewfs()  << std::endl
      << "        max_tau_qpath: " << max_tau_qpath()  << std::endl
      << "        max_tau_spath: " << max_tau_spath()  << std::endl
      << "        max_tau_upath: " << max_tau_upath()  << std::endl
      << "     max_taylor_terms: " << max_taylor_terms()  << std::endl
      << "   max_thermal_coeffs: " << max_thermal_coeffs()  << std::endl
      << "      max_user_levels: " << max_user_levels()  << std::endl
      << "    max_user_obsgeoms: " << max_user_obsgeoms()  << std::endl
      << "     max_user_relazms: " << max_user_relazms()  << std::endl
      << "     max_user_streams: " << max_user_streams()  << std::endl
      << "         maxbandtotal: " << maxbandtotal()  << std::endl
      << "             maxbeams: " << maxbeams()  << std::endl
      << "          maxbrdf_idx: " << maxbrdf_idx()  << std::endl
      << "        maxfinelayers: " << maxfinelayers()  << std::endl
      << "           maxfourier: " << maxfourier()  << std::endl
      << "            maxlayers: " << maxlayers()  << std::endl
      << "           maxmoments: " << maxmoments()  << std::endl
      << "     maxmoments_input: " << maxmoments_input()  << std::endl
      << "       maxsthalf_brdf: " << maxsthalf_brdf()  << std::endl
      << "           maxstreams: " << maxstreams()  << std::endl
      << "         maxstreams_2: " << maxstreams_2()  << std::endl
      << "      maxstreams_brdf: " << maxstreams_brdf()  << std::endl
      << "        maxstreams_p1: " << maxstreams_p1()  << std::endl
      << "   maxstreams_scaling: " << maxstreams_scaling()  << std::endl
      << "             maxtotal: " << maxtotal()  << std::endl
      << "            minus_one: " << minus_one()  << std::endl
      << "            minus_two: " << minus_two()  << std::endl
      << "       modfresnel_idx: " << modfresnel_idx()  << std::endl
      << "       newcmglint_idx: " << newcmglint_idx()  << std::endl
      << "       omega_smallnum: " << omega_smallnum()  << std::endl
      << "                  one: " << one()  << std::endl
      << "                onep5: " << onep5()  << std::endl
      << "                  pi2: " << pi2()  << std::endl
      << "                  pi4: " << pi4()  << std::endl
      << "                  pie: " << pie()  << std::endl
      << "                 pio2: " << pio2()  << std::endl
      << "                 pio4: " << pio4()  << std::endl
      << "              quarter: " << quarter()  << std::endl
      << "           rahman_idx: " << rahman_idx()  << std::endl
      << "        rossthick_idx: " << rossthick_idx()  << std::endl
      << "         rossthin_idx: " << rossthin_idx()  << std::endl
      << "          roujean_idx: " << roujean_idx()  << std::endl
      << "       rtkhotspot_idx: " << rtkhotspot_idx()  << std::endl
      << "             smallnum: " << smallnum()  << std::endl
      << "         snowbrdf_idx: " << snowbrdf_idx()  << std::endl
      << "         taylor_small: " << taylor_small()  << std::endl
      << "                three: " << three()  << std::endl
      << "                  two: " << two()  << std::endl
      << "                upidx: " << upidx()  << std::endl
      << "                 zero: " << zero()  << std::endl;

  }

  Lidort_Pars() { 
    set_lidort_pars(&transfer_struct_c);
  }

private:
  struct lidort_pars_m transfer_struct_c;
};


// Links to type: "brdf_linsup_inputs" from module: "brdf_lin_sup_inputs_def_m" in file: "brdf_lin_sup_inputs_def.f90"
extern "C" {
  void brdf_linsup_inputs_c_alloc_init(struct brdf_linsup_inputs *transfer_struct_c, void **fortran_type_c);
  void brdf_linsup_inputs_c_init_only(struct brdf_linsup_inputs *transfer_struct_c, void **fortran_type_c);
  void brdf_linsup_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_linsup_inputs_c_destroy(void **fortran_type_c);
  
}

struct brdf_linsup_inputs {
  int* bs_do_bsavalue_wf_;
  int bs_do_bsavalue_wf__f_byte_size;

  int* bs_do_kernel_factor_wfs_;
  int bs_do_kernel_factor_wfs__f_shapes[1];
  int bs_do_kernel_factor_wfs__f_byte_size;

  int* bs_do_kernel_params_wfs_;
  int bs_do_kernel_params_wfs__f_shapes[2];
  int bs_do_kernel_params_wfs__f_byte_size;

  int* bs_do_kparams_derivs_;
  int bs_do_kparams_derivs__f_shapes[1];
  int bs_do_kparams_derivs__f_byte_size;

  int* bs_do_windspeed_wf_;
  int bs_do_windspeed_wf__f_byte_size;

  int* bs_do_wsavalue_wf_;
  int bs_do_wsavalue_wf__f_byte_size;

  int* bs_n_kernel_factor_wfs_;
  int bs_n_kernel_factor_wfs__f_byte_size;

  int* bs_n_kernel_params_wfs_;
  int bs_n_kernel_params_wfs__f_byte_size;

  int* bs_n_surface_wfs_;
  int bs_n_surface_wfs__f_byte_size;

  
};

// Links to type: "brdf_linsup_inputs" from module: "brdf_lin_sup_inputs_def_m" in file: "brdf_lin_sup_inputs_def.f90"
class Brdf_Linsup_Inputs : public Brdf_Linsup_Inputs_Base {
public:
  // Allocating constructor
  Brdf_Linsup_Inputs() : Brdf_Linsup_Inputs_Base() {
    brdf_linsup_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Linsup_Inputs(void* allocated_f_type_c) : Brdf_Linsup_Inputs_Base(allocated_f_type_c) {
    brdf_linsup_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Linsup_Inputs() {
    if (owns_pointer)
      brdf_linsup_inputs_c_destroy(&fortran_type_c);
  }

  const bool bs_do_bsavalue_wf() const {
    return *transfer_struct_c.bs_do_bsavalue_wf_ != 0;
  }

  void bs_do_bsavalue_wf(const bool& bs_do_bsavalue_wf_in) {
    *transfer_struct_c.bs_do_bsavalue_wf_ = bs_do_bsavalue_wf_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const blitz::Array<bool, 1> bs_do_kernel_factor_wfs() const {
    blitz::Array<bool,1> as_bool(bs_do_kernel_factor_wfs_.shape());
    as_bool = blitz::where(bs_do_kernel_factor_wfs_ != 0, true, false);
    return as_bool;
  }

  void bs_do_kernel_factor_wfs(const blitz::Array<bool, 1>& bs_do_kernel_factor_wfs_in) {
    blitz::Array<int,1> as_int(bs_do_kernel_factor_wfs_.shape());
    as_int = blitz::where(bs_do_kernel_factor_wfs_in == true, FORTRAN_TRUE_INT, 0);
    bs_do_kernel_factor_wfs_ = as_int;
  }

  
  const blitz::Array<bool, 2> bs_do_kernel_params_wfs() const {
    blitz::Array<bool,2> as_bool(bs_do_kernel_params_wfs_.shape());
    as_bool = blitz::where(bs_do_kernel_params_wfs_ != 0, true, false);
    return as_bool;
  }

  void bs_do_kernel_params_wfs(const blitz::Array<bool, 2>& bs_do_kernel_params_wfs_in) {
    blitz::Array<int,2> as_int(bs_do_kernel_params_wfs_.shape());
    as_int = blitz::where(bs_do_kernel_params_wfs_in == true, FORTRAN_TRUE_INT, 0);
    bs_do_kernel_params_wfs_ = as_int;
  }

  
  const blitz::Array<bool, 1> bs_do_kparams_derivs() const {
    blitz::Array<bool,1> as_bool(bs_do_kparams_derivs_.shape());
    as_bool = blitz::where(bs_do_kparams_derivs_ != 0, true, false);
    return as_bool;
  }

  void bs_do_kparams_derivs(const blitz::Array<bool, 1>& bs_do_kparams_derivs_in) {
    blitz::Array<int,1> as_int(bs_do_kparams_derivs_.shape());
    as_int = blitz::where(bs_do_kparams_derivs_in == true, FORTRAN_TRUE_INT, 0);
    bs_do_kparams_derivs_ = as_int;
  }

  
  const bool bs_do_windspeed_wf() const {
    return *transfer_struct_c.bs_do_windspeed_wf_ != 0;
  }

  void bs_do_windspeed_wf(const bool& bs_do_windspeed_wf_in) {
    *transfer_struct_c.bs_do_windspeed_wf_ = bs_do_windspeed_wf_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_wsavalue_wf() const {
    return *transfer_struct_c.bs_do_wsavalue_wf_ != 0;
  }

  void bs_do_wsavalue_wf(const bool& bs_do_wsavalue_wf_in) {
    *transfer_struct_c.bs_do_wsavalue_wf_ = bs_do_wsavalue_wf_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const int& bs_n_kernel_factor_wfs() const {
    return *transfer_struct_c.bs_n_kernel_factor_wfs_;
  }

  void bs_n_kernel_factor_wfs(const int& bs_n_kernel_factor_wfs_in) {
    *transfer_struct_c.bs_n_kernel_factor_wfs_ = bs_n_kernel_factor_wfs_in;
  }

  
  const int& bs_n_kernel_params_wfs() const {
    return *transfer_struct_c.bs_n_kernel_params_wfs_;
  }

  void bs_n_kernel_params_wfs(const int& bs_n_kernel_params_wfs_in) {
    *transfer_struct_c.bs_n_kernel_params_wfs_ = bs_n_kernel_params_wfs_in;
  }

  
  const int& bs_n_surface_wfs() const {
    return *transfer_struct_c.bs_n_surface_wfs_;
  }

  void bs_n_surface_wfs(const int& bs_n_surface_wfs_in) {
    *transfer_struct_c.bs_n_surface_wfs_ = bs_n_surface_wfs_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Linsup_Inputs:" << std::endl
      << "      bs_do_bsavalue_wf: " << bs_do_bsavalue_wf()  << std::endl
      << "bs_do_kernel_factor_wfs: " << std::endl << bs_do_kernel_factor_wfs()  << std::endl
      << "bs_do_kernel_params_wfs: " << std::endl << bs_do_kernel_params_wfs()  << std::endl
      << "   bs_do_kparams_derivs: " << std::endl << bs_do_kparams_derivs()  << std::endl
      << "     bs_do_windspeed_wf: " << bs_do_windspeed_wf()  << std::endl
      << "      bs_do_wsavalue_wf: " << bs_do_wsavalue_wf()  << std::endl
      << " bs_n_kernel_factor_wfs: " << bs_n_kernel_factor_wfs()  << std::endl
      << " bs_n_kernel_params_wfs: " << bs_n_kernel_params_wfs()  << std::endl
      << "       bs_n_surface_wfs: " << bs_n_surface_wfs()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_do_bsavalue_wf_",sizeof(*transfer_struct_c.bs_do_bsavalue_wf_),transfer_struct_c.bs_do_bsavalue_wf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_kernel_factor_wfs_",sizeof(*transfer_struct_c.bs_do_kernel_factor_wfs_),transfer_struct_c.bs_do_kernel_factor_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_kernel_params_wfs_",sizeof(*transfer_struct_c.bs_do_kernel_params_wfs_),transfer_struct_c.bs_do_kernel_params_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_kparams_derivs_",sizeof(*transfer_struct_c.bs_do_kparams_derivs_),transfer_struct_c.bs_do_kparams_derivs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_windspeed_wf_",sizeof(*transfer_struct_c.bs_do_windspeed_wf_),transfer_struct_c.bs_do_windspeed_wf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_wsavalue_wf_",sizeof(*transfer_struct_c.bs_do_wsavalue_wf_),transfer_struct_c.bs_do_wsavalue_wf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_kernel_factor_wfs_",sizeof(*transfer_struct_c.bs_n_kernel_factor_wfs_),transfer_struct_c.bs_n_kernel_factor_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_kernel_params_wfs_",sizeof(*transfer_struct_c.bs_n_kernel_params_wfs_),transfer_struct_c.bs_n_kernel_params_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_surface_wfs_",sizeof(*transfer_struct_c.bs_n_surface_wfs_),transfer_struct_c.bs_n_surface_wfs__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    bs_do_kernel_factor_wfs_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_do_kernel_factor_wfs_,
      blitz::shape(transfer_struct_c.bs_do_kernel_factor_wfs__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_do_kernel_params_wfs_.reference(blitz::Array<int, 2>(transfer_struct_c.bs_do_kernel_params_wfs_,
      blitz::shape(transfer_struct_c.bs_do_kernel_params_wfs__f_shapes[0],
                   transfer_struct_c.bs_do_kernel_params_wfs__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    bs_do_kparams_derivs_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_do_kparams_derivs_,
      blitz::shape(transfer_struct_c.bs_do_kparams_derivs__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct brdf_linsup_inputs transfer_struct_c;

  blitz::Array<int, 1> bs_do_kernel_factor_wfs_;
  blitz::Array<int, 2> bs_do_kernel_params_wfs_;
  blitz::Array<int, 1> bs_do_kparams_derivs_;
  
};

// Links to type: "brdf_linsup_outputs" from module: "brdf_lin_sup_outputs_def_m" in file: "brdf_lin_sup_outputs_def.f90"
extern "C" {
  void brdf_linsup_outputs_c_alloc_init(struct brdf_linsup_outputs *transfer_struct_c, void **fortran_type_c);
  void brdf_linsup_outputs_c_init_only(struct brdf_linsup_outputs *transfer_struct_c, void **fortran_type_c);
  void brdf_linsup_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_linsup_outputs_c_destroy(void **fortran_type_c);
  
}

struct brdf_linsup_outputs {
  double* bs_ls_brdf_f_;
  int bs_ls_brdf_f__f_shapes[4];
  int bs_ls_brdf_f__f_byte_size;

  double* bs_ls_brdf_f_0_;
  int bs_ls_brdf_f_0__f_shapes[4];
  int bs_ls_brdf_f_0__f_byte_size;

  double* bs_ls_dbounce_brdfunc_;
  int bs_ls_dbounce_brdfunc__f_shapes[4];
  int bs_ls_dbounce_brdfunc__f_byte_size;

  double* bs_ls_emissivity_;
  int bs_ls_emissivity__f_shapes[2];
  int bs_ls_emissivity__f_byte_size;

  double* bs_ls_user_brdf_f_;
  int bs_ls_user_brdf_f__f_shapes[4];
  int bs_ls_user_brdf_f__f_byte_size;

  double* bs_ls_user_brdf_f_0_;
  int bs_ls_user_brdf_f_0__f_shapes[4];
  int bs_ls_user_brdf_f_0__f_byte_size;

  double* bs_ls_user_emissivity_;
  int bs_ls_user_emissivity__f_shapes[2];
  int bs_ls_user_emissivity__f_byte_size;

  
};

// Links to type: "brdf_linsup_outputs" from module: "brdf_lin_sup_outputs_def_m" in file: "brdf_lin_sup_outputs_def.f90"
class Brdf_Linsup_Outputs : public Brdf_Linsup_Outputs_Base {
public:
  // Allocating constructor
  Brdf_Linsup_Outputs() : Brdf_Linsup_Outputs_Base() {
    brdf_linsup_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Linsup_Outputs(void* allocated_f_type_c) : Brdf_Linsup_Outputs_Base(allocated_f_type_c) {
    brdf_linsup_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Linsup_Outputs() {
    if (owns_pointer)
      brdf_linsup_outputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& bs_ls_brdf_f() const {
    return bs_ls_brdf_f_;
  }

  void bs_ls_brdf_f(const blitz::Array<double, 4>& bs_ls_brdf_f_in) {
    bs_ls_brdf_f_ = bs_ls_brdf_f_in;
  }

  
  const blitz::Array<double, 4>& bs_ls_brdf_f_0() const {
    return bs_ls_brdf_f_0_;
  }

  void bs_ls_brdf_f_0(const blitz::Array<double, 4>& bs_ls_brdf_f_0_in) {
    bs_ls_brdf_f_0_ = bs_ls_brdf_f_0_in;
  }

  
  const blitz::Array<double, 4>& bs_ls_dbounce_brdfunc() const {
    return bs_ls_dbounce_brdfunc_;
  }

  void bs_ls_dbounce_brdfunc(const blitz::Array<double, 4>& bs_ls_dbounce_brdfunc_in) {
    bs_ls_dbounce_brdfunc_ = bs_ls_dbounce_brdfunc_in;
  }

  
  const blitz::Array<double, 2>& bs_ls_emissivity() const {
    return bs_ls_emissivity_;
  }

  void bs_ls_emissivity(const blitz::Array<double, 2>& bs_ls_emissivity_in) {
    bs_ls_emissivity_ = bs_ls_emissivity_in;
  }

  
  const blitz::Array<double, 4>& bs_ls_user_brdf_f() const {
    return bs_ls_user_brdf_f_;
  }

  void bs_ls_user_brdf_f(const blitz::Array<double, 4>& bs_ls_user_brdf_f_in) {
    bs_ls_user_brdf_f_ = bs_ls_user_brdf_f_in;
  }

  
  const blitz::Array<double, 4>& bs_ls_user_brdf_f_0() const {
    return bs_ls_user_brdf_f_0_;
  }

  void bs_ls_user_brdf_f_0(const blitz::Array<double, 4>& bs_ls_user_brdf_f_0_in) {
    bs_ls_user_brdf_f_0_ = bs_ls_user_brdf_f_0_in;
  }

  
  const blitz::Array<double, 2>& bs_ls_user_emissivity() const {
    return bs_ls_user_emissivity_;
  }

  void bs_ls_user_emissivity(const blitz::Array<double, 2>& bs_ls_user_emissivity_in) {
    bs_ls_user_emissivity_ = bs_ls_user_emissivity_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Linsup_Outputs:" << std::endl
      << "         bs_ls_brdf_f: " << std::endl << bs_ls_brdf_f()  << std::endl
      << "       bs_ls_brdf_f_0: " << std::endl << bs_ls_brdf_f_0()  << std::endl
      << "bs_ls_dbounce_brdfunc: " << std::endl << bs_ls_dbounce_brdfunc()  << std::endl
      << "     bs_ls_emissivity: " << std::endl << bs_ls_emissivity()  << std::endl
      << "    bs_ls_user_brdf_f: " << std::endl << bs_ls_user_brdf_f()  << std::endl
      << "  bs_ls_user_brdf_f_0: " << std::endl << bs_ls_user_brdf_f_0()  << std::endl
      << "bs_ls_user_emissivity: " << std::endl << bs_ls_user_emissivity()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_ls_brdf_f_",sizeof(*transfer_struct_c.bs_ls_brdf_f_),transfer_struct_c.bs_ls_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_brdf_f_0_",sizeof(*transfer_struct_c.bs_ls_brdf_f_0_),transfer_struct_c.bs_ls_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_dbounce_brdfunc_",sizeof(*transfer_struct_c.bs_ls_dbounce_brdfunc_),transfer_struct_c.bs_ls_dbounce_brdfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_emissivity_",sizeof(*transfer_struct_c.bs_ls_emissivity_),transfer_struct_c.bs_ls_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_user_brdf_f_",sizeof(*transfer_struct_c.bs_ls_user_brdf_f_),transfer_struct_c.bs_ls_user_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_user_brdf_f_0_",sizeof(*transfer_struct_c.bs_ls_user_brdf_f_0_),transfer_struct_c.bs_ls_user_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_ls_user_emissivity_",sizeof(*transfer_struct_c.bs_ls_user_emissivity_),transfer_struct_c.bs_ls_user_emissivity__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    bs_ls_brdf_f_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_brdf_f_,
      blitz::shape(transfer_struct_c.bs_ls_brdf_f__f_shapes[0],
                   transfer_struct_c.bs_ls_brdf_f__f_shapes[1],
                   transfer_struct_c.bs_ls_brdf_f__f_shapes[2],
                   transfer_struct_c.bs_ls_brdf_f__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_brdf_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_brdf_f_0_,
      blitz::shape(transfer_struct_c.bs_ls_brdf_f_0__f_shapes[0],
                   transfer_struct_c.bs_ls_brdf_f_0__f_shapes[1],
                   transfer_struct_c.bs_ls_brdf_f_0__f_shapes[2],
                   transfer_struct_c.bs_ls_brdf_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_dbounce_brdfunc_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_dbounce_brdfunc_,
      blitz::shape(transfer_struct_c.bs_ls_dbounce_brdfunc__f_shapes[0],
                   transfer_struct_c.bs_ls_dbounce_brdfunc__f_shapes[1],
                   transfer_struct_c.bs_ls_dbounce_brdfunc__f_shapes[2],
                   transfer_struct_c.bs_ls_dbounce_brdfunc__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_emissivity_.reference(blitz::Array<double, 2>(transfer_struct_c.bs_ls_emissivity_,
      blitz::shape(transfer_struct_c.bs_ls_emissivity__f_shapes[0],
                   transfer_struct_c.bs_ls_emissivity__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    bs_ls_user_brdf_f_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_user_brdf_f_,
      blitz::shape(transfer_struct_c.bs_ls_user_brdf_f__f_shapes[0],
                   transfer_struct_c.bs_ls_user_brdf_f__f_shapes[1],
                   transfer_struct_c.bs_ls_user_brdf_f__f_shapes[2],
                   transfer_struct_c.bs_ls_user_brdf_f__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_user_brdf_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.bs_ls_user_brdf_f_0_,
      blitz::shape(transfer_struct_c.bs_ls_user_brdf_f_0__f_shapes[0],
                   transfer_struct_c.bs_ls_user_brdf_f_0__f_shapes[1],
                   transfer_struct_c.bs_ls_user_brdf_f_0__f_shapes[2],
                   transfer_struct_c.bs_ls_user_brdf_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    bs_ls_user_emissivity_.reference(blitz::Array<double, 2>(transfer_struct_c.bs_ls_user_emissivity_,
      blitz::shape(transfer_struct_c.bs_ls_user_emissivity__f_shapes[0],
                   transfer_struct_c.bs_ls_user_emissivity__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    
  }

  void link_nested_types() {
    
  }

  struct brdf_linsup_outputs transfer_struct_c;

  blitz::Array<double, 4> bs_ls_brdf_f_;
  blitz::Array<double, 4> bs_ls_brdf_f_0_;
  blitz::Array<double, 4> bs_ls_dbounce_brdfunc_;
  blitz::Array<double, 2> bs_ls_emissivity_;
  blitz::Array<double, 4> bs_ls_user_brdf_f_;
  blitz::Array<double, 4> bs_ls_user_brdf_f_0_;
  blitz::Array<double, 2> bs_ls_user_emissivity_;
  
};

// Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def_m" in file: "brdf_sup_inputs_def.f90"
extern "C" {
  void brdf_sup_inputs_c_alloc_init(struct brdf_sup_inputs *transfer_struct_c, void **fortran_type_c);
  void brdf_sup_inputs_c_init_only(struct brdf_sup_inputs *transfer_struct_c, void **fortran_type_c);
  void brdf_sup_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_sup_inputs_c_destroy(void **fortran_type_c);
  void brdf_sup_inputs_bs_brdf_names_get(void **fortran_type_c, const int* bs_brdf_names_in_shape_1, const int* bs_brdf_names_in_len, const char* bs_brdf_names_in);
  
}

struct brdf_sup_inputs {
  double* bs_beam_szas_;
  int bs_beam_szas__f_shapes[1];
  int bs_beam_szas__f_byte_size;

  double* bs_brdf_factors_;
  int bs_brdf_factors__f_shapes[1];
  int bs_brdf_factors__f_byte_size;

  
  int bs_brdf_names__f_shapes[1];
  int bs_brdf_names__f_len;

  double* bs_brdf_parameters_;
  int bs_brdf_parameters__f_shapes[2];
  int bs_brdf_parameters__f_byte_size;

  double* bs_bsa_value_;
  int bs_bsa_value__f_byte_size;

  int* bs_do_brdf_surface_;
  int bs_do_brdf_surface__f_byte_size;

  int* bs_do_bsa_scaling_;
  int bs_do_bsa_scaling__f_byte_size;

  int* bs_do_directbounce_only_;
  int bs_do_directbounce_only__f_byte_size;

  int* bs_do_doublet_geometry_;
  int bs_do_doublet_geometry__f_byte_size;

  int* bs_do_facetisotropy_;
  int bs_do_facetisotropy__f_byte_size;

  int* bs_do_foamoption_;
  int bs_do_foamoption__f_byte_size;

  int* bs_do_glintshadow_;
  int bs_do_glintshadow__f_byte_size;

  int* bs_do_glitter_msrcorr_;
  int bs_do_glitter_msrcorr__f_byte_size;

  int* bs_do_glitter_msrcorr_dbonly_;
  int bs_do_glitter_msrcorr_dbonly__f_byte_size;

  int* bs_do_newcmglint_;
  int bs_do_newcmglint__f_byte_size;

  int* bs_do_shadow_effect_;
  int bs_do_shadow_effect__f_byte_size;

  int* bs_do_solar_sources_;
  int bs_do_solar_sources__f_byte_size;

  int* bs_do_surface_emission_;
  int bs_do_surface_emission__f_byte_size;

  int* bs_do_user_obsgeoms_;
  int bs_do_user_obsgeoms__f_byte_size;

  int* bs_do_user_streams_;
  int bs_do_user_streams__f_byte_size;

  int* bs_do_wsa_scaling_;
  int bs_do_wsa_scaling__f_byte_size;

  int* bs_do_wsabsa_output_;
  int bs_do_wsabsa_output__f_byte_size;

  int* bs_glitter_msrcorr_nmuquad_;
  int bs_glitter_msrcorr_nmuquad__f_byte_size;

  int* bs_glitter_msrcorr_nphiquad_;
  int bs_glitter_msrcorr_nphiquad__f_byte_size;

  int* bs_glitter_msrcorr_order_;
  int bs_glitter_msrcorr_order__f_byte_size;

  int* bs_lambertian_kernel_flag_;
  int bs_lambertian_kernel_flag__f_shapes[1];
  int bs_lambertian_kernel_flag__f_byte_size;

  int* bs_n_brdf_kernels_;
  int bs_n_brdf_kernels__f_byte_size;

  int* bs_n_brdf_parameters_;
  int bs_n_brdf_parameters__f_shapes[1];
  int bs_n_brdf_parameters__f_byte_size;

  int* bs_n_user_doublets_;
  int bs_n_user_doublets__f_byte_size;

  int* bs_n_user_obsgeoms_;
  int bs_n_user_obsgeoms__f_byte_size;

  int* bs_n_user_relazms_;
  int bs_n_user_relazms__f_byte_size;

  int* bs_n_user_streams_;
  int bs_n_user_streams__f_byte_size;

  int* bs_nbeams_;
  int bs_nbeams__f_byte_size;

  int* bs_nstreams_;
  int bs_nstreams__f_byte_size;

  int* bs_nstreams_brdf_;
  int bs_nstreams_brdf__f_byte_size;

  double* bs_salinity_;
  int bs_salinity__f_byte_size;

  double* bs_user_angles_input_;
  int bs_user_angles_input__f_shapes[1];
  int bs_user_angles_input__f_byte_size;

  double* bs_user_doublets_;
  int bs_user_doublets__f_shapes[2];
  int bs_user_doublets__f_byte_size;

  double* bs_user_obsgeoms_;
  int bs_user_obsgeoms__f_shapes[2];
  int bs_user_obsgeoms__f_byte_size;

  double* bs_user_relazms_;
  int bs_user_relazms__f_shapes[1];
  int bs_user_relazms__f_byte_size;

  double* bs_wavelength_;
  int bs_wavelength__f_byte_size;

  int* bs_which_brdf_;
  int bs_which_brdf__f_shapes[1];
  int bs_which_brdf__f_byte_size;

  double* bs_winddir_;
  int bs_winddir__f_shapes[1];
  int bs_winddir__f_byte_size;

  double* bs_windspeed_;
  int bs_windspeed__f_byte_size;

  double* bs_wsa_value_;
  int bs_wsa_value__f_byte_size;

  
};

// Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def_m" in file: "brdf_sup_inputs_def.f90"
class Brdf_Sup_Inputs : public Brdf_Sup_Inputs_Base {
public:
  // Allocating constructor
  Brdf_Sup_Inputs() : Brdf_Sup_Inputs_Base() {
    brdf_sup_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Sup_Inputs(void* allocated_f_type_c) : Brdf_Sup_Inputs_Base(allocated_f_type_c) {
    brdf_sup_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Sup_Inputs() {
    if (owns_pointer)
      brdf_sup_inputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 1>& bs_beam_szas() const {
    return bs_beam_szas_;
  }

  void bs_beam_szas(const blitz::Array<double, 1>& bs_beam_szas_in) {
    bs_beam_szas_ = bs_beam_szas_in;
  }

  
  const blitz::Array<double, 1>& bs_brdf_factors() const {
    return bs_brdf_factors_;
  }

  void bs_brdf_factors(const blitz::Array<double, 1>& bs_brdf_factors_in) {
    bs_brdf_factors_ = bs_brdf_factors_in;
  }

  
  const std::vector< std::string > bs_brdf_names() const {
    std::vector< std::string > bs_brdf_names_ret;
    blitz::Array<char, 2> bs_brdf_names_lcl = blitz::Array<char, 2>(transfer_struct_c.bs_brdf_names__f_shapes[0], transfer_struct_c.bs_brdf_names__f_len+1, blitz::ColumnMajorArray<2>());
    brdf_sup_inputs_bs_brdf_names_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.bs_brdf_names__f_shapes[0], &transfer_struct_c.bs_brdf_names__f_len, bs_brdf_names_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < bs_brdf_names_lcl.extent(0); dim_0_idx++)
      bs_brdf_names_ret.push_back( std::string(std::string(bs_brdf_names_lcl(dim_0_idx, blitz::Range::all()).begin(), bs_brdf_names_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return bs_brdf_names_ret;
  }

  
  const blitz::Array<double, 2>& bs_brdf_parameters() const {
    return bs_brdf_parameters_;
  }

  void bs_brdf_parameters(const blitz::Array<double, 2>& bs_brdf_parameters_in) {
    bs_brdf_parameters_ = bs_brdf_parameters_in;
  }

  
  const double& bs_bsa_value() const {
    return *transfer_struct_c.bs_bsa_value_;
  }

  void bs_bsa_value(const double& bs_bsa_value_in) {
    *transfer_struct_c.bs_bsa_value_ = bs_bsa_value_in;
  }

  
  const bool bs_do_brdf_surface() const {
    return *transfer_struct_c.bs_do_brdf_surface_ != 0;
  }

  void bs_do_brdf_surface(const bool& bs_do_brdf_surface_in) {
    *transfer_struct_c.bs_do_brdf_surface_ = bs_do_brdf_surface_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_bsa_scaling() const {
    return *transfer_struct_c.bs_do_bsa_scaling_ != 0;
  }

  void bs_do_bsa_scaling(const bool& bs_do_bsa_scaling_in) {
    *transfer_struct_c.bs_do_bsa_scaling_ = bs_do_bsa_scaling_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_directbounce_only() const {
    return *transfer_struct_c.bs_do_directbounce_only_ != 0;
  }

  void bs_do_directbounce_only(const bool& bs_do_directbounce_only_in) {
    *transfer_struct_c.bs_do_directbounce_only_ = bs_do_directbounce_only_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_doublet_geometry() const {
    return *transfer_struct_c.bs_do_doublet_geometry_ != 0;
  }

  void bs_do_doublet_geometry(const bool& bs_do_doublet_geometry_in) {
    *transfer_struct_c.bs_do_doublet_geometry_ = bs_do_doublet_geometry_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_facetisotropy() const {
    return *transfer_struct_c.bs_do_facetisotropy_ != 0;
  }

  void bs_do_facetisotropy(const bool& bs_do_facetisotropy_in) {
    *transfer_struct_c.bs_do_facetisotropy_ = bs_do_facetisotropy_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_foamoption() const {
    return *transfer_struct_c.bs_do_foamoption_ != 0;
  }

  void bs_do_foamoption(const bool& bs_do_foamoption_in) {
    *transfer_struct_c.bs_do_foamoption_ = bs_do_foamoption_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_glintshadow() const {
    return *transfer_struct_c.bs_do_glintshadow_ != 0;
  }

  void bs_do_glintshadow(const bool& bs_do_glintshadow_in) {
    *transfer_struct_c.bs_do_glintshadow_ = bs_do_glintshadow_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_glitter_msrcorr() const {
    return *transfer_struct_c.bs_do_glitter_msrcorr_ != 0;
  }

  void bs_do_glitter_msrcorr(const bool& bs_do_glitter_msrcorr_in) {
    *transfer_struct_c.bs_do_glitter_msrcorr_ = bs_do_glitter_msrcorr_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_glitter_msrcorr_dbonly() const {
    return *transfer_struct_c.bs_do_glitter_msrcorr_dbonly_ != 0;
  }

  void bs_do_glitter_msrcorr_dbonly(const bool& bs_do_glitter_msrcorr_dbonly_in) {
    *transfer_struct_c.bs_do_glitter_msrcorr_dbonly_ = bs_do_glitter_msrcorr_dbonly_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_newcmglint() const {
    return *transfer_struct_c.bs_do_newcmglint_ != 0;
  }

  void bs_do_newcmglint(const bool& bs_do_newcmglint_in) {
    *transfer_struct_c.bs_do_newcmglint_ = bs_do_newcmglint_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_shadow_effect() const {
    return *transfer_struct_c.bs_do_shadow_effect_ != 0;
  }

  void bs_do_shadow_effect(const bool& bs_do_shadow_effect_in) {
    *transfer_struct_c.bs_do_shadow_effect_ = bs_do_shadow_effect_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_solar_sources() const {
    return *transfer_struct_c.bs_do_solar_sources_ != 0;
  }

  void bs_do_solar_sources(const bool& bs_do_solar_sources_in) {
    *transfer_struct_c.bs_do_solar_sources_ = bs_do_solar_sources_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_surface_emission() const {
    return *transfer_struct_c.bs_do_surface_emission_ != 0;
  }

  void bs_do_surface_emission(const bool& bs_do_surface_emission_in) {
    *transfer_struct_c.bs_do_surface_emission_ = bs_do_surface_emission_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_user_obsgeoms() const {
    return *transfer_struct_c.bs_do_user_obsgeoms_ != 0;
  }

  void bs_do_user_obsgeoms(const bool& bs_do_user_obsgeoms_in) {
    *transfer_struct_c.bs_do_user_obsgeoms_ = bs_do_user_obsgeoms_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_user_streams() const {
    return *transfer_struct_c.bs_do_user_streams_ != 0;
  }

  void bs_do_user_streams(const bool& bs_do_user_streams_in) {
    *transfer_struct_c.bs_do_user_streams_ = bs_do_user_streams_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_wsa_scaling() const {
    return *transfer_struct_c.bs_do_wsa_scaling_ != 0;
  }

  void bs_do_wsa_scaling(const bool& bs_do_wsa_scaling_in) {
    *transfer_struct_c.bs_do_wsa_scaling_ = bs_do_wsa_scaling_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool bs_do_wsabsa_output() const {
    return *transfer_struct_c.bs_do_wsabsa_output_ != 0;
  }

  void bs_do_wsabsa_output(const bool& bs_do_wsabsa_output_in) {
    *transfer_struct_c.bs_do_wsabsa_output_ = bs_do_wsabsa_output_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const int& bs_glitter_msrcorr_nmuquad() const {
    return *transfer_struct_c.bs_glitter_msrcorr_nmuquad_;
  }

  void bs_glitter_msrcorr_nmuquad(const int& bs_glitter_msrcorr_nmuquad_in) {
    *transfer_struct_c.bs_glitter_msrcorr_nmuquad_ = bs_glitter_msrcorr_nmuquad_in;
  }

  
  const int& bs_glitter_msrcorr_nphiquad() const {
    return *transfer_struct_c.bs_glitter_msrcorr_nphiquad_;
  }

  void bs_glitter_msrcorr_nphiquad(const int& bs_glitter_msrcorr_nphiquad_in) {
    *transfer_struct_c.bs_glitter_msrcorr_nphiquad_ = bs_glitter_msrcorr_nphiquad_in;
  }

  
  const int& bs_glitter_msrcorr_order() const {
    return *transfer_struct_c.bs_glitter_msrcorr_order_;
  }

  void bs_glitter_msrcorr_order(const int& bs_glitter_msrcorr_order_in) {
    *transfer_struct_c.bs_glitter_msrcorr_order_ = bs_glitter_msrcorr_order_in;
  }

  
  const blitz::Array<bool, 1> bs_lambertian_kernel_flag() const {
    blitz::Array<bool,1> as_bool(bs_lambertian_kernel_flag_.shape());
    as_bool = blitz::where(bs_lambertian_kernel_flag_ != 0, true, false);
    return as_bool;
  }

  void bs_lambertian_kernel_flag(const blitz::Array<bool, 1>& bs_lambertian_kernel_flag_in) {
    blitz::Array<int,1> as_int(bs_lambertian_kernel_flag_.shape());
    as_int = blitz::where(bs_lambertian_kernel_flag_in == true, FORTRAN_TRUE_INT, 0);
    bs_lambertian_kernel_flag_ = as_int;
  }

  
  const int& bs_n_brdf_kernels() const {
    return *transfer_struct_c.bs_n_brdf_kernels_;
  }

  void bs_n_brdf_kernels(const int& bs_n_brdf_kernels_in) {
    *transfer_struct_c.bs_n_brdf_kernels_ = bs_n_brdf_kernels_in;
  }

  
  const blitz::Array<int, 1>& bs_n_brdf_parameters() const {
    return bs_n_brdf_parameters_;
  }

  void bs_n_brdf_parameters(const blitz::Array<int, 1>& bs_n_brdf_parameters_in) {
    bs_n_brdf_parameters_ = bs_n_brdf_parameters_in;
  }

  
  const int& bs_n_user_doublets() const {
    return *transfer_struct_c.bs_n_user_doublets_;
  }

  void bs_n_user_doublets(const int& bs_n_user_doublets_in) {
    *transfer_struct_c.bs_n_user_doublets_ = bs_n_user_doublets_in;
  }

  
  const int& bs_n_user_obsgeoms() const {
    return *transfer_struct_c.bs_n_user_obsgeoms_;
  }

  void bs_n_user_obsgeoms(const int& bs_n_user_obsgeoms_in) {
    *transfer_struct_c.bs_n_user_obsgeoms_ = bs_n_user_obsgeoms_in;
  }

  
  const int& bs_n_user_relazms() const {
    return *transfer_struct_c.bs_n_user_relazms_;
  }

  void bs_n_user_relazms(const int& bs_n_user_relazms_in) {
    *transfer_struct_c.bs_n_user_relazms_ = bs_n_user_relazms_in;
  }

  
  const int& bs_n_user_streams() const {
    return *transfer_struct_c.bs_n_user_streams_;
  }

  void bs_n_user_streams(const int& bs_n_user_streams_in) {
    *transfer_struct_c.bs_n_user_streams_ = bs_n_user_streams_in;
  }

  
  const int& bs_nbeams() const {
    return *transfer_struct_c.bs_nbeams_;
  }

  void bs_nbeams(const int& bs_nbeams_in) {
    *transfer_struct_c.bs_nbeams_ = bs_nbeams_in;
  }

  
  const int& bs_nstreams() const {
    return *transfer_struct_c.bs_nstreams_;
  }

  void bs_nstreams(const int& bs_nstreams_in) {
    *transfer_struct_c.bs_nstreams_ = bs_nstreams_in;
  }

  
  const int& bs_nstreams_brdf() const {
    return *transfer_struct_c.bs_nstreams_brdf_;
  }

  void bs_nstreams_brdf(const int& bs_nstreams_brdf_in) {
    *transfer_struct_c.bs_nstreams_brdf_ = bs_nstreams_brdf_in;
  }

  
  const double& bs_salinity() const {
    return *transfer_struct_c.bs_salinity_;
  }

  void bs_salinity(const double& bs_salinity_in) {
    *transfer_struct_c.bs_salinity_ = bs_salinity_in;
  }

  
  const blitz::Array<double, 1>& bs_user_angles_input() const {
    return bs_user_angles_input_;
  }

  void bs_user_angles_input(const blitz::Array<double, 1>& bs_user_angles_input_in) {
    bs_user_angles_input_ = bs_user_angles_input_in;
  }

  
  const blitz::Array<double, 2>& bs_user_doublets() const {
    return bs_user_doublets_;
  }

  void bs_user_doublets(const blitz::Array<double, 2>& bs_user_doublets_in) {
    bs_user_doublets_ = bs_user_doublets_in;
  }

  
  const blitz::Array<double, 2>& bs_user_obsgeoms() const {
    return bs_user_obsgeoms_;
  }

  void bs_user_obsgeoms(const blitz::Array<double, 2>& bs_user_obsgeoms_in) {
    bs_user_obsgeoms_ = bs_user_obsgeoms_in;
  }

  
  const blitz::Array<double, 1>& bs_user_relazms() const {
    return bs_user_relazms_;
  }

  void bs_user_relazms(const blitz::Array<double, 1>& bs_user_relazms_in) {
    bs_user_relazms_ = bs_user_relazms_in;
  }

  
  const double& bs_wavelength() const {
    return *transfer_struct_c.bs_wavelength_;
  }

  void bs_wavelength(const double& bs_wavelength_in) {
    *transfer_struct_c.bs_wavelength_ = bs_wavelength_in;
  }

  
  const blitz::Array<int, 1>& bs_which_brdf() const {
    return bs_which_brdf_;
  }

  void bs_which_brdf(const blitz::Array<int, 1>& bs_which_brdf_in) {
    bs_which_brdf_ = bs_which_brdf_in;
  }

  
  const blitz::Array<double, 1>& bs_winddir() const {
    return bs_winddir_;
  }

  void bs_winddir(const blitz::Array<double, 1>& bs_winddir_in) {
    bs_winddir_ = bs_winddir_in;
  }

  
  const double& bs_windspeed() const {
    return *transfer_struct_c.bs_windspeed_;
  }

  void bs_windspeed(const double& bs_windspeed_in) {
    *transfer_struct_c.bs_windspeed_ = bs_windspeed_in;
  }

  
  const double& bs_wsa_value() const {
    return *transfer_struct_c.bs_wsa_value_;
  }

  void bs_wsa_value(const double& bs_wsa_value_in) {
    *transfer_struct_c.bs_wsa_value_ = bs_wsa_value_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Sup_Inputs:" << std::endl
      << "                bs_beam_szas: " << std::endl << bs_beam_szas()  << std::endl
      << "             bs_brdf_factors: " << std::endl << bs_brdf_factors()  << std::endl
      << "               bs_brdf_names: " << std::endl;
    std::vector< std::string > bs_brdf_names_lcl = bs_brdf_names();
    for(unsigned int idx = 0; idx < bs_brdf_names_lcl.size(); idx++)
      if ( bs_brdf_names_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << bs_brdf_names_lcl[idx] << "\"" << std::endl;
    output_stream
      << "          bs_brdf_parameters: " << std::endl << bs_brdf_parameters()  << std::endl
      << "                bs_bsa_value: " << bs_bsa_value()  << std::endl
      << "          bs_do_brdf_surface: " << bs_do_brdf_surface()  << std::endl
      << "           bs_do_bsa_scaling: " << bs_do_bsa_scaling()  << std::endl
      << "     bs_do_directbounce_only: " << bs_do_directbounce_only()  << std::endl
      << "      bs_do_doublet_geometry: " << bs_do_doublet_geometry()  << std::endl
      << "         bs_do_facetisotropy: " << bs_do_facetisotropy()  << std::endl
      << "            bs_do_foamoption: " << bs_do_foamoption()  << std::endl
      << "           bs_do_glintshadow: " << bs_do_glintshadow()  << std::endl
      << "       bs_do_glitter_msrcorr: " << bs_do_glitter_msrcorr()  << std::endl
      << "bs_do_glitter_msrcorr_dbonly: " << bs_do_glitter_msrcorr_dbonly()  << std::endl
      << "            bs_do_newcmglint: " << bs_do_newcmglint()  << std::endl
      << "         bs_do_shadow_effect: " << bs_do_shadow_effect()  << std::endl
      << "         bs_do_solar_sources: " << bs_do_solar_sources()  << std::endl
      << "      bs_do_surface_emission: " << bs_do_surface_emission()  << std::endl
      << "         bs_do_user_obsgeoms: " << bs_do_user_obsgeoms()  << std::endl
      << "          bs_do_user_streams: " << bs_do_user_streams()  << std::endl
      << "           bs_do_wsa_scaling: " << bs_do_wsa_scaling()  << std::endl
      << "         bs_do_wsabsa_output: " << bs_do_wsabsa_output()  << std::endl
      << "  bs_glitter_msrcorr_nmuquad: " << bs_glitter_msrcorr_nmuquad()  << std::endl
      << " bs_glitter_msrcorr_nphiquad: " << bs_glitter_msrcorr_nphiquad()  << std::endl
      << "    bs_glitter_msrcorr_order: " << bs_glitter_msrcorr_order()  << std::endl
      << "   bs_lambertian_kernel_flag: " << std::endl << bs_lambertian_kernel_flag()  << std::endl
      << "           bs_n_brdf_kernels: " << bs_n_brdf_kernels()  << std::endl
      << "        bs_n_brdf_parameters: " << std::endl << bs_n_brdf_parameters()  << std::endl
      << "          bs_n_user_doublets: " << bs_n_user_doublets()  << std::endl
      << "          bs_n_user_obsgeoms: " << bs_n_user_obsgeoms()  << std::endl
      << "           bs_n_user_relazms: " << bs_n_user_relazms()  << std::endl
      << "           bs_n_user_streams: " << bs_n_user_streams()  << std::endl
      << "                   bs_nbeams: " << bs_nbeams()  << std::endl
      << "                 bs_nstreams: " << bs_nstreams()  << std::endl
      << "            bs_nstreams_brdf: " << bs_nstreams_brdf()  << std::endl
      << "                 bs_salinity: " << bs_salinity()  << std::endl
      << "        bs_user_angles_input: " << std::endl << bs_user_angles_input()  << std::endl
      << "            bs_user_doublets: " << std::endl << bs_user_doublets()  << std::endl
      << "            bs_user_obsgeoms: " << std::endl << bs_user_obsgeoms()  << std::endl
      << "             bs_user_relazms: " << std::endl << bs_user_relazms()  << std::endl
      << "               bs_wavelength: " << bs_wavelength()  << std::endl
      << "               bs_which_brdf: " << std::endl << bs_which_brdf()  << std::endl
      << "                  bs_winddir: " << std::endl << bs_winddir()  << std::endl
      << "                bs_windspeed: " << bs_windspeed()  << std::endl
      << "                bs_wsa_value: " << bs_wsa_value()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_beam_szas_",sizeof(*transfer_struct_c.bs_beam_szas_),transfer_struct_c.bs_beam_szas__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_brdf_factors_",sizeof(*transfer_struct_c.bs_brdf_factors_),transfer_struct_c.bs_brdf_factors__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_brdf_parameters_",sizeof(*transfer_struct_c.bs_brdf_parameters_),transfer_struct_c.bs_brdf_parameters__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_bsa_value_",sizeof(*transfer_struct_c.bs_bsa_value_),transfer_struct_c.bs_bsa_value__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_brdf_surface_",sizeof(*transfer_struct_c.bs_do_brdf_surface_),transfer_struct_c.bs_do_brdf_surface__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_bsa_scaling_",sizeof(*transfer_struct_c.bs_do_bsa_scaling_),transfer_struct_c.bs_do_bsa_scaling__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_directbounce_only_",sizeof(*transfer_struct_c.bs_do_directbounce_only_),transfer_struct_c.bs_do_directbounce_only__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_doublet_geometry_",sizeof(*transfer_struct_c.bs_do_doublet_geometry_),transfer_struct_c.bs_do_doublet_geometry__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_facetisotropy_",sizeof(*transfer_struct_c.bs_do_facetisotropy_),transfer_struct_c.bs_do_facetisotropy__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_foamoption_",sizeof(*transfer_struct_c.bs_do_foamoption_),transfer_struct_c.bs_do_foamoption__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_glintshadow_",sizeof(*transfer_struct_c.bs_do_glintshadow_),transfer_struct_c.bs_do_glintshadow__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_glitter_msrcorr_",sizeof(*transfer_struct_c.bs_do_glitter_msrcorr_),transfer_struct_c.bs_do_glitter_msrcorr__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_glitter_msrcorr_dbonly_",sizeof(*transfer_struct_c.bs_do_glitter_msrcorr_dbonly_),transfer_struct_c.bs_do_glitter_msrcorr_dbonly__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_newcmglint_",sizeof(*transfer_struct_c.bs_do_newcmglint_),transfer_struct_c.bs_do_newcmglint__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_shadow_effect_",sizeof(*transfer_struct_c.bs_do_shadow_effect_),transfer_struct_c.bs_do_shadow_effect__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_solar_sources_",sizeof(*transfer_struct_c.bs_do_solar_sources_),transfer_struct_c.bs_do_solar_sources__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_surface_emission_",sizeof(*transfer_struct_c.bs_do_surface_emission_),transfer_struct_c.bs_do_surface_emission__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_user_obsgeoms_",sizeof(*transfer_struct_c.bs_do_user_obsgeoms_),transfer_struct_c.bs_do_user_obsgeoms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_user_streams_",sizeof(*transfer_struct_c.bs_do_user_streams_),transfer_struct_c.bs_do_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_wsa_scaling_",sizeof(*transfer_struct_c.bs_do_wsa_scaling_),transfer_struct_c.bs_do_wsa_scaling__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_do_wsabsa_output_",sizeof(*transfer_struct_c.bs_do_wsabsa_output_),transfer_struct_c.bs_do_wsabsa_output__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_glitter_msrcorr_nmuquad_",sizeof(*transfer_struct_c.bs_glitter_msrcorr_nmuquad_),transfer_struct_c.bs_glitter_msrcorr_nmuquad__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_glitter_msrcorr_nphiquad_",sizeof(*transfer_struct_c.bs_glitter_msrcorr_nphiquad_),transfer_struct_c.bs_glitter_msrcorr_nphiquad__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_glitter_msrcorr_order_",sizeof(*transfer_struct_c.bs_glitter_msrcorr_order_),transfer_struct_c.bs_glitter_msrcorr_order__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_lambertian_kernel_flag_",sizeof(*transfer_struct_c.bs_lambertian_kernel_flag_),transfer_struct_c.bs_lambertian_kernel_flag__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_brdf_kernels_",sizeof(*transfer_struct_c.bs_n_brdf_kernels_),transfer_struct_c.bs_n_brdf_kernels__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_brdf_parameters_",sizeof(*transfer_struct_c.bs_n_brdf_parameters_),transfer_struct_c.bs_n_brdf_parameters__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_user_doublets_",sizeof(*transfer_struct_c.bs_n_user_doublets_),transfer_struct_c.bs_n_user_doublets__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_user_obsgeoms_",sizeof(*transfer_struct_c.bs_n_user_obsgeoms_),transfer_struct_c.bs_n_user_obsgeoms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_user_relazms_",sizeof(*transfer_struct_c.bs_n_user_relazms_),transfer_struct_c.bs_n_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_n_user_streams_",sizeof(*transfer_struct_c.bs_n_user_streams_),transfer_struct_c.bs_n_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_nbeams_",sizeof(*transfer_struct_c.bs_nbeams_),transfer_struct_c.bs_nbeams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_nstreams_",sizeof(*transfer_struct_c.bs_nstreams_),transfer_struct_c.bs_nstreams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_nstreams_brdf_",sizeof(*transfer_struct_c.bs_nstreams_brdf_),transfer_struct_c.bs_nstreams_brdf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_salinity_",sizeof(*transfer_struct_c.bs_salinity_),transfer_struct_c.bs_salinity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_angles_input_",sizeof(*transfer_struct_c.bs_user_angles_input_),transfer_struct_c.bs_user_angles_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_doublets_",sizeof(*transfer_struct_c.bs_user_doublets_),transfer_struct_c.bs_user_doublets__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_obsgeoms_",sizeof(*transfer_struct_c.bs_user_obsgeoms_),transfer_struct_c.bs_user_obsgeoms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_relazms_",sizeof(*transfer_struct_c.bs_user_relazms_),transfer_struct_c.bs_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_wavelength_",sizeof(*transfer_struct_c.bs_wavelength_),transfer_struct_c.bs_wavelength__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_which_brdf_",sizeof(*transfer_struct_c.bs_which_brdf_),transfer_struct_c.bs_which_brdf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_winddir_",sizeof(*transfer_struct_c.bs_winddir_),transfer_struct_c.bs_winddir__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_windspeed_",sizeof(*transfer_struct_c.bs_windspeed_),transfer_struct_c.bs_windspeed__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_wsa_value_",sizeof(*transfer_struct_c.bs_wsa_value_),transfer_struct_c.bs_wsa_value__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    bs_beam_szas_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_beam_szas_,
      blitz::shape(transfer_struct_c.bs_beam_szas__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_brdf_factors_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_brdf_factors_,
      blitz::shape(transfer_struct_c.bs_brdf_factors__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_brdf_parameters_.reference(blitz::Array<double, 2>(transfer_struct_c.bs_brdf_parameters_,
      blitz::shape(transfer_struct_c.bs_brdf_parameters__f_shapes[0],
                   transfer_struct_c.bs_brdf_parameters__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    bs_lambertian_kernel_flag_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_lambertian_kernel_flag_,
      blitz::shape(transfer_struct_c.bs_lambertian_kernel_flag__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_n_brdf_parameters_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_n_brdf_parameters_,
      blitz::shape(transfer_struct_c.bs_n_brdf_parameters__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_user_angles_input_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_user_angles_input_,
      blitz::shape(transfer_struct_c.bs_user_angles_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_user_doublets_.reference(blitz::Array<double, 2>(transfer_struct_c.bs_user_doublets_,
      blitz::shape(transfer_struct_c.bs_user_doublets__f_shapes[0],
                   transfer_struct_c.bs_user_doublets__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    bs_user_obsgeoms_.reference(blitz::Array<double, 2>(transfer_struct_c.bs_user_obsgeoms_,
      blitz::shape(transfer_struct_c.bs_user_obsgeoms__f_shapes[0],
                   transfer_struct_c.bs_user_obsgeoms__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    bs_user_relazms_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_user_relazms_,
      blitz::shape(transfer_struct_c.bs_user_relazms__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_which_brdf_.reference(blitz::Array<int, 1>(transfer_struct_c.bs_which_brdf_,
      blitz::shape(transfer_struct_c.bs_which_brdf__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_winddir_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_winddir_,
      blitz::shape(transfer_struct_c.bs_winddir__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct brdf_sup_inputs transfer_struct_c;

  blitz::Array<double, 1> bs_beam_szas_;
  blitz::Array<double, 1> bs_brdf_factors_;
  blitz::Array<double, 2> bs_brdf_parameters_;
  blitz::Array<int, 1> bs_lambertian_kernel_flag_;
  blitz::Array<int, 1> bs_n_brdf_parameters_;
  blitz::Array<double, 1> bs_user_angles_input_;
  blitz::Array<double, 2> bs_user_doublets_;
  blitz::Array<double, 2> bs_user_obsgeoms_;
  blitz::Array<double, 1> bs_user_relazms_;
  blitz::Array<int, 1> bs_which_brdf_;
  blitz::Array<double, 1> bs_winddir_;
  
};

// Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.f90"
extern "C" {
  void brdf_sup_outputs_c_alloc_init(struct brdf_sup_outputs *transfer_struct_c, void **fortran_type_c);
  void brdf_sup_outputs_c_init_only(struct brdf_sup_outputs *transfer_struct_c, void **fortran_type_c);
  void brdf_sup_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_sup_outputs_c_destroy(void **fortran_type_c);
  
}

struct brdf_sup_outputs {
  double* bs_brdf_f_;
  int bs_brdf_f__f_shapes[3];
  int bs_brdf_f__f_byte_size;

  double* bs_brdf_f_0_;
  int bs_brdf_f_0__f_shapes[3];
  int bs_brdf_f_0__f_byte_size;

  double* bs_bsa_calculated_;
  int bs_bsa_calculated__f_byte_size;

  double* bs_bsa_kernels_;
  int bs_bsa_kernels__f_shapes[1];
  int bs_bsa_kernels__f_byte_size;

  double* bs_dbounce_brdfunc_;
  int bs_dbounce_brdfunc__f_shapes[3];
  int bs_dbounce_brdfunc__f_byte_size;

  double* bs_emissivity_;
  int bs_emissivity__f_shapes[1];
  int bs_emissivity__f_byte_size;

  double* bs_user_brdf_f_;
  int bs_user_brdf_f__f_shapes[3];
  int bs_user_brdf_f__f_byte_size;

  double* bs_user_brdf_f_0_;
  int bs_user_brdf_f_0__f_shapes[3];
  int bs_user_brdf_f_0__f_byte_size;

  double* bs_user_emissivity_;
  int bs_user_emissivity__f_shapes[1];
  int bs_user_emissivity__f_byte_size;

  double* bs_wsa_calculated_;
  int bs_wsa_calculated__f_byte_size;

  double* bs_wsa_kernels_;
  int bs_wsa_kernels__f_shapes[1];
  int bs_wsa_kernels__f_byte_size;

  
};

// Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.f90"
class Brdf_Sup_Outputs : public Brdf_Sup_Outputs_Base {
public:
  // Allocating constructor
  Brdf_Sup_Outputs() : Brdf_Sup_Outputs_Base() {
    brdf_sup_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Sup_Outputs(void* allocated_f_type_c) : Brdf_Sup_Outputs_Base(allocated_f_type_c) {
    brdf_sup_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Sup_Outputs() {
    if (owns_pointer)
      brdf_sup_outputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& bs_brdf_f() const {
    return bs_brdf_f_;
  }

  void bs_brdf_f(const blitz::Array<double, 3>& bs_brdf_f_in) {
    bs_brdf_f_ = bs_brdf_f_in;
  }

  
  const blitz::Array<double, 3>& bs_brdf_f_0() const {
    return bs_brdf_f_0_;
  }

  void bs_brdf_f_0(const blitz::Array<double, 3>& bs_brdf_f_0_in) {
    bs_brdf_f_0_ = bs_brdf_f_0_in;
  }

  
  const double& bs_bsa_calculated() const {
    return *transfer_struct_c.bs_bsa_calculated_;
  }

  void bs_bsa_calculated(const double& bs_bsa_calculated_in) {
    *transfer_struct_c.bs_bsa_calculated_ = bs_bsa_calculated_in;
  }

  
  const blitz::Array<double, 1>& bs_bsa_kernels() const {
    return bs_bsa_kernels_;
  }

  void bs_bsa_kernels(const blitz::Array<double, 1>& bs_bsa_kernels_in) {
    bs_bsa_kernels_ = bs_bsa_kernels_in;
  }

  
  const blitz::Array<double, 3>& bs_dbounce_brdfunc() const {
    return bs_dbounce_brdfunc_;
  }

  void bs_dbounce_brdfunc(const blitz::Array<double, 3>& bs_dbounce_brdfunc_in) {
    bs_dbounce_brdfunc_ = bs_dbounce_brdfunc_in;
  }

  
  const blitz::Array<double, 1>& bs_emissivity() const {
    return bs_emissivity_;
  }

  void bs_emissivity(const blitz::Array<double, 1>& bs_emissivity_in) {
    bs_emissivity_ = bs_emissivity_in;
  }

  
  const blitz::Array<double, 3>& bs_user_brdf_f() const {
    return bs_user_brdf_f_;
  }

  void bs_user_brdf_f(const blitz::Array<double, 3>& bs_user_brdf_f_in) {
    bs_user_brdf_f_ = bs_user_brdf_f_in;
  }

  
  const blitz::Array<double, 3>& bs_user_brdf_f_0() const {
    return bs_user_brdf_f_0_;
  }

  void bs_user_brdf_f_0(const blitz::Array<double, 3>& bs_user_brdf_f_0_in) {
    bs_user_brdf_f_0_ = bs_user_brdf_f_0_in;
  }

  
  const blitz::Array<double, 1>& bs_user_emissivity() const {
    return bs_user_emissivity_;
  }

  void bs_user_emissivity(const blitz::Array<double, 1>& bs_user_emissivity_in) {
    bs_user_emissivity_ = bs_user_emissivity_in;
  }

  
  const double& bs_wsa_calculated() const {
    return *transfer_struct_c.bs_wsa_calculated_;
  }

  void bs_wsa_calculated(const double& bs_wsa_calculated_in) {
    *transfer_struct_c.bs_wsa_calculated_ = bs_wsa_calculated_in;
  }

  
  const blitz::Array<double, 1>& bs_wsa_kernels() const {
    return bs_wsa_kernels_;
  }

  void bs_wsa_kernels(const blitz::Array<double, 1>& bs_wsa_kernels_in) {
    bs_wsa_kernels_ = bs_wsa_kernels_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Sup_Outputs:" << std::endl
      << "         bs_brdf_f: " << std::endl << bs_brdf_f()  << std::endl
      << "       bs_brdf_f_0: " << std::endl << bs_brdf_f_0()  << std::endl
      << " bs_bsa_calculated: " << bs_bsa_calculated()  << std::endl
      << "    bs_bsa_kernels: " << std::endl << bs_bsa_kernels()  << std::endl
      << "bs_dbounce_brdfunc: " << std::endl << bs_dbounce_brdfunc()  << std::endl
      << "     bs_emissivity: " << std::endl << bs_emissivity()  << std::endl
      << "    bs_user_brdf_f: " << std::endl << bs_user_brdf_f()  << std::endl
      << "  bs_user_brdf_f_0: " << std::endl << bs_user_brdf_f_0()  << std::endl
      << "bs_user_emissivity: " << std::endl << bs_user_emissivity()  << std::endl
      << " bs_wsa_calculated: " << bs_wsa_calculated()  << std::endl
      << "    bs_wsa_kernels: " << std::endl << bs_wsa_kernels()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_brdf_f_",sizeof(*transfer_struct_c.bs_brdf_f_),transfer_struct_c.bs_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_brdf_f_0_",sizeof(*transfer_struct_c.bs_brdf_f_0_),transfer_struct_c.bs_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_bsa_calculated_",sizeof(*transfer_struct_c.bs_bsa_calculated_),transfer_struct_c.bs_bsa_calculated__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_bsa_kernels_",sizeof(*transfer_struct_c.bs_bsa_kernels_),transfer_struct_c.bs_bsa_kernels__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_dbounce_brdfunc_",sizeof(*transfer_struct_c.bs_dbounce_brdfunc_),transfer_struct_c.bs_dbounce_brdfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_emissivity_",sizeof(*transfer_struct_c.bs_emissivity_),transfer_struct_c.bs_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_brdf_f_",sizeof(*transfer_struct_c.bs_user_brdf_f_),transfer_struct_c.bs_user_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_brdf_f_0_",sizeof(*transfer_struct_c.bs_user_brdf_f_0_),transfer_struct_c.bs_user_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_user_emissivity_",sizeof(*transfer_struct_c.bs_user_emissivity_),transfer_struct_c.bs_user_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_wsa_calculated_",sizeof(*transfer_struct_c.bs_wsa_calculated_),transfer_struct_c.bs_wsa_calculated__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_wsa_kernels_",sizeof(*transfer_struct_c.bs_wsa_kernels_),transfer_struct_c.bs_wsa_kernels__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    bs_brdf_f_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_brdf_f_,
      blitz::shape(transfer_struct_c.bs_brdf_f__f_shapes[0],
                   transfer_struct_c.bs_brdf_f__f_shapes[1],
                   transfer_struct_c.bs_brdf_f__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_brdf_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_brdf_f_0_,
      blitz::shape(transfer_struct_c.bs_brdf_f_0__f_shapes[0],
                   transfer_struct_c.bs_brdf_f_0__f_shapes[1],
                   transfer_struct_c.bs_brdf_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_bsa_kernels_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_bsa_kernels_,
      blitz::shape(transfer_struct_c.bs_bsa_kernels__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_dbounce_brdfunc_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_dbounce_brdfunc_,
      blitz::shape(transfer_struct_c.bs_dbounce_brdfunc__f_shapes[0],
                   transfer_struct_c.bs_dbounce_brdfunc__f_shapes[1],
                   transfer_struct_c.bs_dbounce_brdfunc__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_emissivity_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_emissivity_,
      blitz::shape(transfer_struct_c.bs_emissivity__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_user_brdf_f_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_user_brdf_f_,
      blitz::shape(transfer_struct_c.bs_user_brdf_f__f_shapes[0],
                   transfer_struct_c.bs_user_brdf_f__f_shapes[1],
                   transfer_struct_c.bs_user_brdf_f__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_user_brdf_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.bs_user_brdf_f_0_,
      blitz::shape(transfer_struct_c.bs_user_brdf_f_0__f_shapes[0],
                   transfer_struct_c.bs_user_brdf_f_0__f_shapes[1],
                   transfer_struct_c.bs_user_brdf_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    bs_user_emissivity_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_user_emissivity_,
      blitz::shape(transfer_struct_c.bs_user_emissivity__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    bs_wsa_kernels_.reference(blitz::Array<double, 1>(transfer_struct_c.bs_wsa_kernels_,
      blitz::shape(transfer_struct_c.bs_wsa_kernels__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct brdf_sup_outputs transfer_struct_c;

  blitz::Array<double, 3> bs_brdf_f_;
  blitz::Array<double, 3> bs_brdf_f_0_;
  blitz::Array<double, 1> bs_bsa_kernels_;
  blitz::Array<double, 3> bs_dbounce_brdfunc_;
  blitz::Array<double, 1> bs_emissivity_;
  blitz::Array<double, 3> bs_user_brdf_f_;
  blitz::Array<double, 3> bs_user_brdf_f_0_;
  blitz::Array<double, 1> bs_user_emissivity_;
  blitz::Array<double, 1> bs_wsa_kernels_;
  
};

// Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.f90"
extern "C" {
  void brdf_input_exception_handling_c_alloc_init(struct brdf_input_exception_handling *transfer_struct_c, void **fortran_type_c);
  void brdf_input_exception_handling_c_init_only(struct brdf_input_exception_handling *transfer_struct_c, void **fortran_type_c);
  void brdf_input_exception_handling_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_input_exception_handling_c_destroy(void **fortran_type_c);
  void brdf_input_exception_handling_bs_inputactions_get(void **fortran_type_c, const int* bs_inputactions_in_shape_1, const int* bs_inputactions_in_len, const char* bs_inputactions_in);
  void brdf_input_exception_handling_bs_inputmessages_get(void **fortran_type_c, const int* bs_inputmessages_in_shape_1, const int* bs_inputmessages_in_len, const char* bs_inputmessages_in);
  
}

struct brdf_input_exception_handling {
  
  int bs_inputactions__f_shapes[1];
  int bs_inputactions__f_len;

  
  int bs_inputmessages__f_shapes[1];
  int bs_inputmessages__f_len;

  int* bs_ninputmessages_;
  int bs_ninputmessages__f_byte_size;

  int* bs_status_inputread_;
  int bs_status_inputread__f_byte_size;

  
};

// Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.f90"
class Brdf_Input_Exception_Handling : public Brdf_Input_Exception_Handling_Base {
public:
  // Allocating constructor
  Brdf_Input_Exception_Handling() : Brdf_Input_Exception_Handling_Base() {
    brdf_input_exception_handling_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Input_Exception_Handling(void* allocated_f_type_c) : Brdf_Input_Exception_Handling_Base(allocated_f_type_c) {
    brdf_input_exception_handling_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Input_Exception_Handling() {
    if (owns_pointer)
      brdf_input_exception_handling_c_destroy(&fortran_type_c);
  }

  const std::vector< std::string > bs_inputactions() const {
    std::vector< std::string > bs_inputactions_ret;
    blitz::Array<char, 2> bs_inputactions_lcl = blitz::Array<char, 2>(transfer_struct_c.bs_inputactions__f_shapes[0], transfer_struct_c.bs_inputactions__f_len+1, blitz::ColumnMajorArray<2>());
    brdf_input_exception_handling_bs_inputactions_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.bs_inputactions__f_shapes[0], &transfer_struct_c.bs_inputactions__f_len, bs_inputactions_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < bs_inputactions_lcl.extent(0); dim_0_idx++)
      bs_inputactions_ret.push_back( std::string(std::string(bs_inputactions_lcl(dim_0_idx, blitz::Range::all()).begin(), bs_inputactions_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return bs_inputactions_ret;
  }

  
  const std::vector< std::string > bs_inputmessages() const {
    std::vector< std::string > bs_inputmessages_ret;
    blitz::Array<char, 2> bs_inputmessages_lcl = blitz::Array<char, 2>(transfer_struct_c.bs_inputmessages__f_shapes[0], transfer_struct_c.bs_inputmessages__f_len+1, blitz::ColumnMajorArray<2>());
    brdf_input_exception_handling_bs_inputmessages_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.bs_inputmessages__f_shapes[0], &transfer_struct_c.bs_inputmessages__f_len, bs_inputmessages_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < bs_inputmessages_lcl.extent(0); dim_0_idx++)
      bs_inputmessages_ret.push_back( std::string(std::string(bs_inputmessages_lcl(dim_0_idx, blitz::Range::all()).begin(), bs_inputmessages_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return bs_inputmessages_ret;
  }

  
  const int& bs_ninputmessages() const {
    return *transfer_struct_c.bs_ninputmessages_;
  }

  void bs_ninputmessages(const int& bs_ninputmessages_in) {
    *transfer_struct_c.bs_ninputmessages_ = bs_ninputmessages_in;
  }

  
  const int& bs_status_inputread() const {
    return *transfer_struct_c.bs_status_inputread_;
  }

  void bs_status_inputread(const int& bs_status_inputread_in) {
    *transfer_struct_c.bs_status_inputread_ = bs_status_inputread_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Input_Exception_Handling:" << std::endl
      << "    bs_inputactions: " << std::endl;
    std::vector< std::string > bs_inputactions_lcl = bs_inputactions();
    for(unsigned int idx = 0; idx < bs_inputactions_lcl.size(); idx++)
      if ( bs_inputactions_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << bs_inputactions_lcl[idx] << "\"" << std::endl;
    output_stream
      << "   bs_inputmessages: " << std::endl;
    std::vector< std::string > bs_inputmessages_lcl = bs_inputmessages();
    for(unsigned int idx = 0; idx < bs_inputmessages_lcl.size(); idx++)
      if ( bs_inputmessages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << bs_inputmessages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "  bs_ninputmessages: " << bs_ninputmessages()  << std::endl
      << "bs_status_inputread: " << bs_status_inputread()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_ninputmessages_",sizeof(*transfer_struct_c.bs_ninputmessages_),transfer_struct_c.bs_ninputmessages__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_status_inputread_",sizeof(*transfer_struct_c.bs_status_inputread_),transfer_struct_c.bs_status_inputread__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct brdf_input_exception_handling transfer_struct_c;

  
};

// Links to type: "brdf_output_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.f90"
extern "C" {
  void brdf_output_exception_handling_c_alloc_init(struct brdf_output_exception_handling *transfer_struct_c, void **fortran_type_c);
  void brdf_output_exception_handling_c_init_only(struct brdf_output_exception_handling *transfer_struct_c, void **fortran_type_c);
  void brdf_output_exception_handling_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void brdf_output_exception_handling_c_destroy(void **fortran_type_c);
  void brdf_output_exception_handling_bs_outputmessages_get(void **fortran_type_c, const int* bs_outputmessages_in_shape_1, const int* bs_outputmessages_in_len, const char* bs_outputmessages_in);
  
}

struct brdf_output_exception_handling {
  int* bs_noutputmessages_;
  int bs_noutputmessages__f_byte_size;

  
  int bs_outputmessages__f_shapes[1];
  int bs_outputmessages__f_len;

  int* bs_status_output_;
  int bs_status_output__f_byte_size;

  
};

// Links to type: "brdf_output_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.f90"
class Brdf_Output_Exception_Handling : public Brdf_Output_Exception_Handling_Base {
public:
  // Allocating constructor
  Brdf_Output_Exception_Handling() : Brdf_Output_Exception_Handling_Base() {
    brdf_output_exception_handling_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Brdf_Output_Exception_Handling(void* allocated_f_type_c) : Brdf_Output_Exception_Handling_Base(allocated_f_type_c) {
    brdf_output_exception_handling_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Brdf_Output_Exception_Handling() {
    if (owns_pointer)
      brdf_output_exception_handling_c_destroy(&fortran_type_c);
  }

  const int& bs_noutputmessages() const {
    return *transfer_struct_c.bs_noutputmessages_;
  }

  void bs_noutputmessages(const int& bs_noutputmessages_in) {
    *transfer_struct_c.bs_noutputmessages_ = bs_noutputmessages_in;
  }

  
  const std::vector< std::string > bs_outputmessages() const {
    std::vector< std::string > bs_outputmessages_ret;
    blitz::Array<char, 2> bs_outputmessages_lcl = blitz::Array<char, 2>(transfer_struct_c.bs_outputmessages__f_shapes[0], transfer_struct_c.bs_outputmessages__f_len+1, blitz::ColumnMajorArray<2>());
    brdf_output_exception_handling_bs_outputmessages_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.bs_outputmessages__f_shapes[0], &transfer_struct_c.bs_outputmessages__f_len, bs_outputmessages_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < bs_outputmessages_lcl.extent(0); dim_0_idx++)
      bs_outputmessages_ret.push_back( std::string(std::string(bs_outputmessages_lcl(dim_0_idx, blitz::Range::all()).begin(), bs_outputmessages_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return bs_outputmessages_ret;
  }

  
  const int& bs_status_output() const {
    return *transfer_struct_c.bs_status_output_;
  }

  void bs_status_output(const int& bs_status_output_in) {
    *transfer_struct_c.bs_status_output_ = bs_status_output_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Brdf_Output_Exception_Handling:" << std::endl
      << "bs_noutputmessages: " << bs_noutputmessages()  << std::endl
      << " bs_outputmessages: " << std::endl;
    std::vector< std::string > bs_outputmessages_lcl = bs_outputmessages();
    for(unsigned int idx = 0; idx < bs_outputmessages_lcl.size(); idx++)
      if ( bs_outputmessages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << bs_outputmessages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "  bs_status_output: " << bs_status_output()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("bs_noutputmessages_",sizeof(*transfer_struct_c.bs_noutputmessages_),transfer_struct_c.bs_noutputmessages__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("bs_status_output_",sizeof(*transfer_struct_c.bs_status_output_),transfer_struct_c.bs_status_output__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct brdf_output_exception_handling transfer_struct_c;

  
};

// Links to type: "sleave_sup_inputs" from module: "sleave_sup_inputs_def_m" in file: "sleave_sup_inputs_def.f90"
extern "C" {
  void sleave_sup_inputs_c_alloc_init(struct sleave_sup_inputs *transfer_struct_c, void **fortran_type_c);
  void sleave_sup_inputs_c_init_only(struct sleave_sup_inputs *transfer_struct_c, void **fortran_type_c);
  void sleave_sup_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void sleave_sup_inputs_c_destroy(void **fortran_type_c);
  void sleave_sup_inputs_sl_sleave_datapath_get(void **fortran_type_c, const int* sl_sleave_datapath_in_len, const char* sl_sleave_datapath_in);
  
}

struct sleave_sup_inputs {
  int* sl_azimuthdep_;
  int sl_azimuthdep__f_byte_size;

  double* sl_beam_szas_;
  int sl_beam_szas__f_shapes[1];
  int sl_beam_szas__f_byte_size;

  double* sl_chlorconc_;
  int sl_chlorconc__f_byte_size;

  int* sl_do_doublet_geometry_;
  int sl_do_doublet_geometry__f_byte_size;

  int* sl_do_exact_;
  int sl_do_exact__f_byte_size;

  int* sl_do_exactonly_;
  int sl_do_exactonly__f_byte_size;

  int* sl_do_facetisotropy_;
  int sl_do_facetisotropy__f_byte_size;

  int* sl_do_fluorescence_;
  int sl_do_fluorescence__f_byte_size;

  int* sl_do_foamoption_;
  int sl_do_foamoption__f_byte_size;

  int* sl_do_fourier_output_;
  int sl_do_fourier_output__f_byte_size;

  int* sl_do_glintshadow_;
  int sl_do_glintshadow__f_byte_size;

  int* sl_do_isotropic_;
  int sl_do_isotropic__f_byte_size;

  int* sl_do_roughsurface_;
  int sl_do_roughsurface__f_byte_size;

  int* sl_do_sleaving_;
  int sl_do_sleaving__f_byte_size;

  int* sl_do_solar_sources_;
  int sl_do_solar_sources__f_byte_size;

  int* sl_do_user_obsgeoms_;
  int sl_do_user_obsgeoms__f_byte_size;

  int* sl_do_user_streams_;
  int sl_do_user_streams__f_byte_size;

  double* sl_fl_amplitude755_;
  int sl_fl_amplitude755__f_byte_size;

  int* sl_fl_do_datagaussian_;
  int sl_fl_do_datagaussian__f_byte_size;

  int* sl_fl_epoch_;
  int sl_fl_epoch__f_shapes[1];
  int sl_fl_epoch__f_byte_size;

  double* sl_fl_inputgaussians_;
  int sl_fl_inputgaussians__f_shapes[2];
  int sl_fl_inputgaussians__f_byte_size;

  double* sl_fl_latitude_;
  int sl_fl_latitude__f_byte_size;

  double* sl_fl_longitude_;
  int sl_fl_longitude__f_byte_size;

  double* sl_fl_wavelength_;
  int sl_fl_wavelength__f_byte_size;

  int* sl_n_user_doublets_;
  int sl_n_user_doublets__f_byte_size;

  int* sl_n_user_obsgeoms_;
  int sl_n_user_obsgeoms__f_byte_size;

  int* sl_n_user_relazms_;
  int sl_n_user_relazms__f_byte_size;

  int* sl_n_user_streams_;
  int sl_n_user_streams__f_byte_size;

  int* sl_nbeams_;
  int sl_nbeams__f_byte_size;

  int* sl_nstreams_;
  int sl_nstreams__f_byte_size;

  double* sl_salinity_;
  int sl_salinity__f_byte_size;

  
  int sl_sleave_datapath__f_len;

  double* sl_user_angles_input_;
  int sl_user_angles_input__f_shapes[1];
  int sl_user_angles_input__f_byte_size;

  double* sl_user_doublets_;
  int sl_user_doublets__f_shapes[2];
  int sl_user_doublets__f_byte_size;

  double* sl_user_obsgeoms_;
  int sl_user_obsgeoms__f_shapes[2];
  int sl_user_obsgeoms__f_byte_size;

  double* sl_user_relazms_;
  int sl_user_relazms__f_shapes[1];
  int sl_user_relazms__f_byte_size;

  double* sl_wavelength_;
  int sl_wavelength__f_byte_size;

  double* sl_winddir_;
  int sl_winddir__f_shapes[1];
  int sl_winddir__f_byte_size;

  double* sl_windspeed_;
  int sl_windspeed__f_byte_size;

  
};

// Links to type: "sleave_sup_inputs" from module: "sleave_sup_inputs_def_m" in file: "sleave_sup_inputs_def.f90"
class Sleave_Sup_Inputs : public Sleave_Sup_Inputs_Base {
public:
  // Allocating constructor
  Sleave_Sup_Inputs() : Sleave_Sup_Inputs_Base() {
    sleave_sup_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Sleave_Sup_Inputs(void* allocated_f_type_c) : Sleave_Sup_Inputs_Base(allocated_f_type_c) {
    sleave_sup_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Sleave_Sup_Inputs() {
    if (owns_pointer)
      sleave_sup_inputs_c_destroy(&fortran_type_c);
  }

  const bool sl_azimuthdep() const {
    return *transfer_struct_c.sl_azimuthdep_ != 0;
  }

  void sl_azimuthdep(const bool& sl_azimuthdep_in) {
    *transfer_struct_c.sl_azimuthdep_ = sl_azimuthdep_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const blitz::Array<double, 1>& sl_beam_szas() const {
    return sl_beam_szas_;
  }

  void sl_beam_szas(const blitz::Array<double, 1>& sl_beam_szas_in) {
    sl_beam_szas_ = sl_beam_szas_in;
  }

  
  const double& sl_chlorconc() const {
    return *transfer_struct_c.sl_chlorconc_;
  }

  void sl_chlorconc(const double& sl_chlorconc_in) {
    *transfer_struct_c.sl_chlorconc_ = sl_chlorconc_in;
  }

  
  const bool sl_do_doublet_geometry() const {
    return *transfer_struct_c.sl_do_doublet_geometry_ != 0;
  }

  void sl_do_doublet_geometry(const bool& sl_do_doublet_geometry_in) {
    *transfer_struct_c.sl_do_doublet_geometry_ = sl_do_doublet_geometry_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_exact() const {
    return *transfer_struct_c.sl_do_exact_ != 0;
  }

  void sl_do_exact(const bool& sl_do_exact_in) {
    *transfer_struct_c.sl_do_exact_ = sl_do_exact_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_exactonly() const {
    return *transfer_struct_c.sl_do_exactonly_ != 0;
  }

  void sl_do_exactonly(const bool& sl_do_exactonly_in) {
    *transfer_struct_c.sl_do_exactonly_ = sl_do_exactonly_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_facetisotropy() const {
    return *transfer_struct_c.sl_do_facetisotropy_ != 0;
  }

  void sl_do_facetisotropy(const bool& sl_do_facetisotropy_in) {
    *transfer_struct_c.sl_do_facetisotropy_ = sl_do_facetisotropy_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_fluorescence() const {
    return *transfer_struct_c.sl_do_fluorescence_ != 0;
  }

  void sl_do_fluorescence(const bool& sl_do_fluorescence_in) {
    *transfer_struct_c.sl_do_fluorescence_ = sl_do_fluorescence_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_foamoption() const {
    return *transfer_struct_c.sl_do_foamoption_ != 0;
  }

  void sl_do_foamoption(const bool& sl_do_foamoption_in) {
    *transfer_struct_c.sl_do_foamoption_ = sl_do_foamoption_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_fourier_output() const {
    return *transfer_struct_c.sl_do_fourier_output_ != 0;
  }

  void sl_do_fourier_output(const bool& sl_do_fourier_output_in) {
    *transfer_struct_c.sl_do_fourier_output_ = sl_do_fourier_output_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_glintshadow() const {
    return *transfer_struct_c.sl_do_glintshadow_ != 0;
  }

  void sl_do_glintshadow(const bool& sl_do_glintshadow_in) {
    *transfer_struct_c.sl_do_glintshadow_ = sl_do_glintshadow_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_isotropic() const {
    return *transfer_struct_c.sl_do_isotropic_ != 0;
  }

  void sl_do_isotropic(const bool& sl_do_isotropic_in) {
    *transfer_struct_c.sl_do_isotropic_ = sl_do_isotropic_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_roughsurface() const {
    return *transfer_struct_c.sl_do_roughsurface_ != 0;
  }

  void sl_do_roughsurface(const bool& sl_do_roughsurface_in) {
    *transfer_struct_c.sl_do_roughsurface_ = sl_do_roughsurface_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_sleaving() const {
    return *transfer_struct_c.sl_do_sleaving_ != 0;
  }

  void sl_do_sleaving(const bool& sl_do_sleaving_in) {
    *transfer_struct_c.sl_do_sleaving_ = sl_do_sleaving_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_solar_sources() const {
    return *transfer_struct_c.sl_do_solar_sources_ != 0;
  }

  void sl_do_solar_sources(const bool& sl_do_solar_sources_in) {
    *transfer_struct_c.sl_do_solar_sources_ = sl_do_solar_sources_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_user_obsgeoms() const {
    return *transfer_struct_c.sl_do_user_obsgeoms_ != 0;
  }

  void sl_do_user_obsgeoms(const bool& sl_do_user_obsgeoms_in) {
    *transfer_struct_c.sl_do_user_obsgeoms_ = sl_do_user_obsgeoms_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool sl_do_user_streams() const {
    return *transfer_struct_c.sl_do_user_streams_ != 0;
  }

  void sl_do_user_streams(const bool& sl_do_user_streams_in) {
    *transfer_struct_c.sl_do_user_streams_ = sl_do_user_streams_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const double& sl_fl_amplitude755() const {
    return *transfer_struct_c.sl_fl_amplitude755_;
  }

  void sl_fl_amplitude755(const double& sl_fl_amplitude755_in) {
    *transfer_struct_c.sl_fl_amplitude755_ = sl_fl_amplitude755_in;
  }

  
  const bool sl_fl_do_datagaussian() const {
    return *transfer_struct_c.sl_fl_do_datagaussian_ != 0;
  }

  void sl_fl_do_datagaussian(const bool& sl_fl_do_datagaussian_in) {
    *transfer_struct_c.sl_fl_do_datagaussian_ = sl_fl_do_datagaussian_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const blitz::Array<int, 1>& sl_fl_epoch() const {
    return sl_fl_epoch_;
  }

  void sl_fl_epoch(const blitz::Array<int, 1>& sl_fl_epoch_in) {
    sl_fl_epoch_ = sl_fl_epoch_in;
  }

  
  const blitz::Array<double, 2>& sl_fl_inputgaussians() const {
    return sl_fl_inputgaussians_;
  }

  void sl_fl_inputgaussians(const blitz::Array<double, 2>& sl_fl_inputgaussians_in) {
    sl_fl_inputgaussians_ = sl_fl_inputgaussians_in;
  }

  
  const double& sl_fl_latitude() const {
    return *transfer_struct_c.sl_fl_latitude_;
  }

  void sl_fl_latitude(const double& sl_fl_latitude_in) {
    *transfer_struct_c.sl_fl_latitude_ = sl_fl_latitude_in;
  }

  
  const double& sl_fl_longitude() const {
    return *transfer_struct_c.sl_fl_longitude_;
  }

  void sl_fl_longitude(const double& sl_fl_longitude_in) {
    *transfer_struct_c.sl_fl_longitude_ = sl_fl_longitude_in;
  }

  
  const double& sl_fl_wavelength() const {
    return *transfer_struct_c.sl_fl_wavelength_;
  }

  void sl_fl_wavelength(const double& sl_fl_wavelength_in) {
    *transfer_struct_c.sl_fl_wavelength_ = sl_fl_wavelength_in;
  }

  
  const int& sl_n_user_doublets() const {
    return *transfer_struct_c.sl_n_user_doublets_;
  }

  void sl_n_user_doublets(const int& sl_n_user_doublets_in) {
    *transfer_struct_c.sl_n_user_doublets_ = sl_n_user_doublets_in;
  }

  
  const int& sl_n_user_obsgeoms() const {
    return *transfer_struct_c.sl_n_user_obsgeoms_;
  }

  void sl_n_user_obsgeoms(const int& sl_n_user_obsgeoms_in) {
    *transfer_struct_c.sl_n_user_obsgeoms_ = sl_n_user_obsgeoms_in;
  }

  
  const int& sl_n_user_relazms() const {
    return *transfer_struct_c.sl_n_user_relazms_;
  }

  void sl_n_user_relazms(const int& sl_n_user_relazms_in) {
    *transfer_struct_c.sl_n_user_relazms_ = sl_n_user_relazms_in;
  }

  
  const int& sl_n_user_streams() const {
    return *transfer_struct_c.sl_n_user_streams_;
  }

  void sl_n_user_streams(const int& sl_n_user_streams_in) {
    *transfer_struct_c.sl_n_user_streams_ = sl_n_user_streams_in;
  }

  
  const int& sl_nbeams() const {
    return *transfer_struct_c.sl_nbeams_;
  }

  void sl_nbeams(const int& sl_nbeams_in) {
    *transfer_struct_c.sl_nbeams_ = sl_nbeams_in;
  }

  
  const int& sl_nstreams() const {
    return *transfer_struct_c.sl_nstreams_;
  }

  void sl_nstreams(const int& sl_nstreams_in) {
    *transfer_struct_c.sl_nstreams_ = sl_nstreams_in;
  }

  
  const double& sl_salinity() const {
    return *transfer_struct_c.sl_salinity_;
  }

  void sl_salinity(const double& sl_salinity_in) {
    *transfer_struct_c.sl_salinity_ = sl_salinity_in;
  }

  
  const std::string sl_sleave_datapath() const {
    std::string sl_sleave_datapath_ret;
    blitz::Array<char, 1> sl_sleave_datapath_lcl = blitz::Array<char, 1>(transfer_struct_c.sl_sleave_datapath__f_len+1, blitz::ColumnMajorArray<1>());
    sleave_sup_inputs_sl_sleave_datapath_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.sl_sleave_datapath__f_len, sl_sleave_datapath_lcl.dataFirst());
    sl_sleave_datapath_ret = ( std::string(std::string(sl_sleave_datapath_lcl(blitz::Range::all()).begin(), sl_sleave_datapath_lcl(blitz::Range::all()).end()).c_str()) );
    return sl_sleave_datapath_ret;
  }

  
  const blitz::Array<double, 1>& sl_user_angles_input() const {
    return sl_user_angles_input_;
  }

  void sl_user_angles_input(const blitz::Array<double, 1>& sl_user_angles_input_in) {
    sl_user_angles_input_ = sl_user_angles_input_in;
  }

  
  const blitz::Array<double, 2>& sl_user_doublets() const {
    return sl_user_doublets_;
  }

  void sl_user_doublets(const blitz::Array<double, 2>& sl_user_doublets_in) {
    sl_user_doublets_ = sl_user_doublets_in;
  }

  
  const blitz::Array<double, 2>& sl_user_obsgeoms() const {
    return sl_user_obsgeoms_;
  }

  void sl_user_obsgeoms(const blitz::Array<double, 2>& sl_user_obsgeoms_in) {
    sl_user_obsgeoms_ = sl_user_obsgeoms_in;
  }

  
  const blitz::Array<double, 1>& sl_user_relazms() const {
    return sl_user_relazms_;
  }

  void sl_user_relazms(const blitz::Array<double, 1>& sl_user_relazms_in) {
    sl_user_relazms_ = sl_user_relazms_in;
  }

  
  const double& sl_wavelength() const {
    return *transfer_struct_c.sl_wavelength_;
  }

  void sl_wavelength(const double& sl_wavelength_in) {
    *transfer_struct_c.sl_wavelength_ = sl_wavelength_in;
  }

  
  const blitz::Array<double, 1>& sl_winddir() const {
    return sl_winddir_;
  }

  void sl_winddir(const blitz::Array<double, 1>& sl_winddir_in) {
    sl_winddir_ = sl_winddir_in;
  }

  
  const double& sl_windspeed() const {
    return *transfer_struct_c.sl_windspeed_;
  }

  void sl_windspeed(const double& sl_windspeed_in) {
    *transfer_struct_c.sl_windspeed_ = sl_windspeed_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Sleave_Sup_Inputs:" << std::endl
      << "         sl_azimuthdep: " << sl_azimuthdep()  << std::endl
      << "          sl_beam_szas: " << std::endl << sl_beam_szas()  << std::endl
      << "          sl_chlorconc: " << sl_chlorconc()  << std::endl
      << "sl_do_doublet_geometry: " << sl_do_doublet_geometry()  << std::endl
      << "           sl_do_exact: " << sl_do_exact()  << std::endl
      << "       sl_do_exactonly: " << sl_do_exactonly()  << std::endl
      << "   sl_do_facetisotropy: " << sl_do_facetisotropy()  << std::endl
      << "    sl_do_fluorescence: " << sl_do_fluorescence()  << std::endl
      << "      sl_do_foamoption: " << sl_do_foamoption()  << std::endl
      << "  sl_do_fourier_output: " << sl_do_fourier_output()  << std::endl
      << "     sl_do_glintshadow: " << sl_do_glintshadow()  << std::endl
      << "       sl_do_isotropic: " << sl_do_isotropic()  << std::endl
      << "    sl_do_roughsurface: " << sl_do_roughsurface()  << std::endl
      << "        sl_do_sleaving: " << sl_do_sleaving()  << std::endl
      << "   sl_do_solar_sources: " << sl_do_solar_sources()  << std::endl
      << "   sl_do_user_obsgeoms: " << sl_do_user_obsgeoms()  << std::endl
      << "    sl_do_user_streams: " << sl_do_user_streams()  << std::endl
      << "    sl_fl_amplitude755: " << sl_fl_amplitude755()  << std::endl
      << " sl_fl_do_datagaussian: " << sl_fl_do_datagaussian()  << std::endl
      << "           sl_fl_epoch: " << std::endl << sl_fl_epoch()  << std::endl
      << "  sl_fl_inputgaussians: " << std::endl << sl_fl_inputgaussians()  << std::endl
      << "        sl_fl_latitude: " << sl_fl_latitude()  << std::endl
      << "       sl_fl_longitude: " << sl_fl_longitude()  << std::endl
      << "      sl_fl_wavelength: " << sl_fl_wavelength()  << std::endl
      << "    sl_n_user_doublets: " << sl_n_user_doublets()  << std::endl
      << "    sl_n_user_obsgeoms: " << sl_n_user_obsgeoms()  << std::endl
      << "     sl_n_user_relazms: " << sl_n_user_relazms()  << std::endl
      << "     sl_n_user_streams: " << sl_n_user_streams()  << std::endl
      << "             sl_nbeams: " << sl_nbeams()  << std::endl
      << "           sl_nstreams: " << sl_nstreams()  << std::endl
      << "           sl_salinity: " << sl_salinity()  << std::endl
      << "    sl_sleave_datapath: " << "\"" << sl_sleave_datapath() << "\"" << std::endl
      << "  sl_user_angles_input: " << std::endl << sl_user_angles_input()  << std::endl
      << "      sl_user_doublets: " << std::endl << sl_user_doublets()  << std::endl
      << "      sl_user_obsgeoms: " << std::endl << sl_user_obsgeoms()  << std::endl
      << "       sl_user_relazms: " << std::endl << sl_user_relazms()  << std::endl
      << "         sl_wavelength: " << sl_wavelength()  << std::endl
      << "            sl_winddir: " << std::endl << sl_winddir()  << std::endl
      << "          sl_windspeed: " << sl_windspeed()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("sl_azimuthdep_",sizeof(*transfer_struct_c.sl_azimuthdep_),transfer_struct_c.sl_azimuthdep__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_beam_szas_",sizeof(*transfer_struct_c.sl_beam_szas_),transfer_struct_c.sl_beam_szas__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_chlorconc_",sizeof(*transfer_struct_c.sl_chlorconc_),transfer_struct_c.sl_chlorconc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_doublet_geometry_",sizeof(*transfer_struct_c.sl_do_doublet_geometry_),transfer_struct_c.sl_do_doublet_geometry__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_exact_",sizeof(*transfer_struct_c.sl_do_exact_),transfer_struct_c.sl_do_exact__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_exactonly_",sizeof(*transfer_struct_c.sl_do_exactonly_),transfer_struct_c.sl_do_exactonly__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_facetisotropy_",sizeof(*transfer_struct_c.sl_do_facetisotropy_),transfer_struct_c.sl_do_facetisotropy__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_fluorescence_",sizeof(*transfer_struct_c.sl_do_fluorescence_),transfer_struct_c.sl_do_fluorescence__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_foamoption_",sizeof(*transfer_struct_c.sl_do_foamoption_),transfer_struct_c.sl_do_foamoption__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_fourier_output_",sizeof(*transfer_struct_c.sl_do_fourier_output_),transfer_struct_c.sl_do_fourier_output__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_glintshadow_",sizeof(*transfer_struct_c.sl_do_glintshadow_),transfer_struct_c.sl_do_glintshadow__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_isotropic_",sizeof(*transfer_struct_c.sl_do_isotropic_),transfer_struct_c.sl_do_isotropic__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_roughsurface_",sizeof(*transfer_struct_c.sl_do_roughsurface_),transfer_struct_c.sl_do_roughsurface__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_sleaving_",sizeof(*transfer_struct_c.sl_do_sleaving_),transfer_struct_c.sl_do_sleaving__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_solar_sources_",sizeof(*transfer_struct_c.sl_do_solar_sources_),transfer_struct_c.sl_do_solar_sources__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_user_obsgeoms_",sizeof(*transfer_struct_c.sl_do_user_obsgeoms_),transfer_struct_c.sl_do_user_obsgeoms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_do_user_streams_",sizeof(*transfer_struct_c.sl_do_user_streams_),transfer_struct_c.sl_do_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_fl_amplitude755_",sizeof(*transfer_struct_c.sl_fl_amplitude755_),transfer_struct_c.sl_fl_amplitude755__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_fl_do_datagaussian_",sizeof(*transfer_struct_c.sl_fl_do_datagaussian_),transfer_struct_c.sl_fl_do_datagaussian__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_fl_epoch_",sizeof(*transfer_struct_c.sl_fl_epoch_),transfer_struct_c.sl_fl_epoch__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_fl_inputgaussians_",sizeof(*transfer_struct_c.sl_fl_inputgaussians_),transfer_struct_c.sl_fl_inputgaussians__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_fl_latitude_",sizeof(*transfer_struct_c.sl_fl_latitude_),transfer_struct_c.sl_fl_latitude__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_fl_longitude_",sizeof(*transfer_struct_c.sl_fl_longitude_),transfer_struct_c.sl_fl_longitude__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_fl_wavelength_",sizeof(*transfer_struct_c.sl_fl_wavelength_),transfer_struct_c.sl_fl_wavelength__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_n_user_doublets_",sizeof(*transfer_struct_c.sl_n_user_doublets_),transfer_struct_c.sl_n_user_doublets__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_n_user_obsgeoms_",sizeof(*transfer_struct_c.sl_n_user_obsgeoms_),transfer_struct_c.sl_n_user_obsgeoms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_n_user_relazms_",sizeof(*transfer_struct_c.sl_n_user_relazms_),transfer_struct_c.sl_n_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_n_user_streams_",sizeof(*transfer_struct_c.sl_n_user_streams_),transfer_struct_c.sl_n_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_nbeams_",sizeof(*transfer_struct_c.sl_nbeams_),transfer_struct_c.sl_nbeams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_nstreams_",sizeof(*transfer_struct_c.sl_nstreams_),transfer_struct_c.sl_nstreams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_salinity_",sizeof(*transfer_struct_c.sl_salinity_),transfer_struct_c.sl_salinity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_user_angles_input_",sizeof(*transfer_struct_c.sl_user_angles_input_),transfer_struct_c.sl_user_angles_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_user_doublets_",sizeof(*transfer_struct_c.sl_user_doublets_),transfer_struct_c.sl_user_doublets__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_user_obsgeoms_",sizeof(*transfer_struct_c.sl_user_obsgeoms_),transfer_struct_c.sl_user_obsgeoms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_user_relazms_",sizeof(*transfer_struct_c.sl_user_relazms_),transfer_struct_c.sl_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_wavelength_",sizeof(*transfer_struct_c.sl_wavelength_),transfer_struct_c.sl_wavelength__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_winddir_",sizeof(*transfer_struct_c.sl_winddir_),transfer_struct_c.sl_winddir__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("sl_windspeed_",sizeof(*transfer_struct_c.sl_windspeed_),transfer_struct_c.sl_windspeed__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    sl_beam_szas_.reference(blitz::Array<double, 1>(transfer_struct_c.sl_beam_szas_,
      blitz::shape(transfer_struct_c.sl_beam_szas__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    sl_fl_epoch_.reference(blitz::Array<int, 1>(transfer_struct_c.sl_fl_epoch_,
      blitz::shape(transfer_struct_c.sl_fl_epoch__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    sl_fl_inputgaussians_.reference(blitz::Array<double, 2>(transfer_struct_c.sl_fl_inputgaussians_,
      blitz::shape(transfer_struct_c.sl_fl_inputgaussians__f_shapes[0],
                   transfer_struct_c.sl_fl_inputgaussians__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    sl_user_angles_input_.reference(blitz::Array<double, 1>(transfer_struct_c.sl_user_angles_input_,
      blitz::shape(transfer_struct_c.sl_user_angles_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    sl_user_doublets_.reference(blitz::Array<double, 2>(transfer_struct_c.sl_user_doublets_,
      blitz::shape(transfer_struct_c.sl_user_doublets__f_shapes[0],
                   transfer_struct_c.sl_user_doublets__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    sl_user_obsgeoms_.reference(blitz::Array<double, 2>(transfer_struct_c.sl_user_obsgeoms_,
      blitz::shape(transfer_struct_c.sl_user_obsgeoms__f_shapes[0],
                   transfer_struct_c.sl_user_obsgeoms__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    sl_user_relazms_.reference(blitz::Array<double, 1>(transfer_struct_c.sl_user_relazms_,
      blitz::shape(transfer_struct_c.sl_user_relazms__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    sl_winddir_.reference(blitz::Array<double, 1>(transfer_struct_c.sl_winddir_,
      blitz::shape(transfer_struct_c.sl_winddir__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct sleave_sup_inputs transfer_struct_c;

  blitz::Array<double, 1> sl_beam_szas_;
  blitz::Array<int, 1> sl_fl_epoch_;
  blitz::Array<double, 2> sl_fl_inputgaussians_;
  blitz::Array<double, 1> sl_user_angles_input_;
  blitz::Array<double, 2> sl_user_doublets_;
  blitz::Array<double, 2> sl_user_obsgeoms_;
  blitz::Array<double, 1> sl_user_relazms_;
  blitz::Array<double, 1> sl_winddir_;
  
};

// Links to type: "lidort_fixed_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
extern "C" {
  void lidort_fixed_lincontrol_c_alloc_init(struct lidort_fixed_lincontrol *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_lincontrol_c_init_only(struct lidort_fixed_lincontrol *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_lincontrol_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_lincontrol_c_destroy(void **fortran_type_c);
  void fixed_lincontrol_ts_columnwf_names_get(void **fortran_type_c, const int* ts_columnwf_names_in_shape_1, const int* ts_columnwf_names_in_len, const char* ts_columnwf_names_in);
  void fixed_lincontrol_ts_profilewf_names_get(void **fortran_type_c, const int* ts_profilewf_names_in_shape_1, const int* ts_profilewf_names_in_len, const char* ts_profilewf_names_in);
  
}

struct lidort_fixed_lincontrol {
  
  int ts_columnwf_names__f_shapes[1];
  int ts_columnwf_names__f_len;

  int* ts_layer_vary_flag_;
  int ts_layer_vary_flag__f_shapes[1];
  int ts_layer_vary_flag__f_byte_size;

  int* ts_layer_vary_number_;
  int ts_layer_vary_number__f_shapes[1];
  int ts_layer_vary_number__f_byte_size;

  int* ts_n_sleave_wfs_;
  int ts_n_sleave_wfs__f_byte_size;

  int* ts_n_surface_wfs_;
  int ts_n_surface_wfs__f_byte_size;

  int* ts_n_totalcolumn_wfs_;
  int ts_n_totalcolumn_wfs__f_byte_size;

  
  int ts_profilewf_names__f_shapes[1];
  int ts_profilewf_names__f_len;

  
};

// Links to type: "lidort_fixed_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
class Lidort_Fixed_Lincontrol : public Spurr_Fixed_Lincontrol_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Lincontrol() : Spurr_Fixed_Lincontrol_Base() {
    lidort_fixed_lincontrol_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Lincontrol(void* allocated_f_type_c) : Spurr_Fixed_Lincontrol_Base(allocated_f_type_c) {
    lidort_fixed_lincontrol_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Lincontrol() {
    if (owns_pointer)
      lidort_fixed_lincontrol_c_destroy(&fortran_type_c);
  }

  const std::vector< std::string > ts_columnwf_names() const {
    std::vector< std::string > ts_columnwf_names_ret;
    blitz::Array<char, 2> ts_columnwf_names_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_columnwf_names__f_shapes[0], transfer_struct_c.ts_columnwf_names__f_len+1, blitz::ColumnMajorArray<2>());
    fixed_lincontrol_ts_columnwf_names_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_columnwf_names__f_shapes[0], &transfer_struct_c.ts_columnwf_names__f_len, ts_columnwf_names_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_columnwf_names_lcl.extent(0); dim_0_idx++)
      ts_columnwf_names_ret.push_back( std::string(std::string(ts_columnwf_names_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_columnwf_names_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_columnwf_names_ret;
  }

  
  const blitz::Array<bool, 1> ts_layer_vary_flag() const {
    blitz::Array<bool,1> as_bool(ts_layer_vary_flag_.shape());
    as_bool = blitz::where(ts_layer_vary_flag_ != 0, true, false);
    return as_bool;
  }

  void ts_layer_vary_flag(const blitz::Array<bool, 1>& ts_layer_vary_flag_in) {
    blitz::Array<int,1> as_int(ts_layer_vary_flag_.shape());
    as_int = blitz::where(ts_layer_vary_flag_in == true, FORTRAN_TRUE_INT, 0);
    ts_layer_vary_flag_ = as_int;
  }

  
  const blitz::Array<int, 1>& ts_layer_vary_number() const {
    return ts_layer_vary_number_;
  }

  void ts_layer_vary_number(const blitz::Array<int, 1>& ts_layer_vary_number_in) {
    ts_layer_vary_number_ = ts_layer_vary_number_in;
  }

  
  const int& ts_n_sleave_wfs() const {
    return *transfer_struct_c.ts_n_sleave_wfs_;
  }

  void ts_n_sleave_wfs(const int& ts_n_sleave_wfs_in) {
    *transfer_struct_c.ts_n_sleave_wfs_ = ts_n_sleave_wfs_in;
  }

  
  const int& ts_n_surface_wfs() const {
    return *transfer_struct_c.ts_n_surface_wfs_;
  }

  void ts_n_surface_wfs(const int& ts_n_surface_wfs_in) {
    *transfer_struct_c.ts_n_surface_wfs_ = ts_n_surface_wfs_in;
  }

  
  const int& ts_n_totalcolumn_wfs() const {
    return *transfer_struct_c.ts_n_totalcolumn_wfs_;
  }

  void ts_n_totalcolumn_wfs(const int& ts_n_totalcolumn_wfs_in) {
    *transfer_struct_c.ts_n_totalcolumn_wfs_ = ts_n_totalcolumn_wfs_in;
  }

  
  const std::vector< std::string > ts_profilewf_names() const {
    std::vector< std::string > ts_profilewf_names_ret;
    blitz::Array<char, 2> ts_profilewf_names_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_profilewf_names__f_shapes[0], transfer_struct_c.ts_profilewf_names__f_len+1, blitz::ColumnMajorArray<2>());
    fixed_lincontrol_ts_profilewf_names_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_profilewf_names__f_shapes[0], &transfer_struct_c.ts_profilewf_names__f_len, ts_profilewf_names_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_profilewf_names_lcl.extent(0); dim_0_idx++)
      ts_profilewf_names_ret.push_back( std::string(std::string(ts_profilewf_names_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_profilewf_names_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_profilewf_names_ret;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Lincontrol:" << std::endl
      << "   ts_columnwf_names: " << std::endl;
    std::vector< std::string > ts_columnwf_names_lcl = ts_columnwf_names();
    for(unsigned int idx = 0; idx < ts_columnwf_names_lcl.size(); idx++)
      if ( ts_columnwf_names_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_columnwf_names_lcl[idx] << "\"" << std::endl;
    output_stream
      << "  ts_layer_vary_flag: " << std::endl << ts_layer_vary_flag()  << std::endl
      << "ts_layer_vary_number: " << std::endl << ts_layer_vary_number()  << std::endl
      << "     ts_n_sleave_wfs: " << ts_n_sleave_wfs()  << std::endl
      << "    ts_n_surface_wfs: " << ts_n_surface_wfs()  << std::endl
      << "ts_n_totalcolumn_wfs: " << ts_n_totalcolumn_wfs()  << std::endl
      << "  ts_profilewf_names: " << std::endl;
    std::vector< std::string > ts_profilewf_names_lcl = ts_profilewf_names();
    for(unsigned int idx = 0; idx < ts_profilewf_names_lcl.size(); idx++)
      if ( ts_profilewf_names_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_profilewf_names_lcl[idx] << "\"" << std::endl;
  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_layer_vary_flag_",sizeof(*transfer_struct_c.ts_layer_vary_flag_),transfer_struct_c.ts_layer_vary_flag__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_layer_vary_number_",sizeof(*transfer_struct_c.ts_layer_vary_number_),transfer_struct_c.ts_layer_vary_number__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_sleave_wfs_",sizeof(*transfer_struct_c.ts_n_sleave_wfs_),transfer_struct_c.ts_n_sleave_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_surface_wfs_",sizeof(*transfer_struct_c.ts_n_surface_wfs_),transfer_struct_c.ts_n_surface_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_totalcolumn_wfs_",sizeof(*transfer_struct_c.ts_n_totalcolumn_wfs_),transfer_struct_c.ts_n_totalcolumn_wfs__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_layer_vary_flag_.reference(blitz::Array<int, 1>(transfer_struct_c.ts_layer_vary_flag_,
      blitz::shape(transfer_struct_c.ts_layer_vary_flag__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_layer_vary_number_.reference(blitz::Array<int, 1>(transfer_struct_c.ts_layer_vary_number_,
      blitz::shape(transfer_struct_c.ts_layer_vary_number__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_lincontrol transfer_struct_c;

  blitz::Array<int, 1> ts_layer_vary_flag_;
  blitz::Array<int, 1> ts_layer_vary_number_;
  
};

// Links to type: "lidort_fixed_linoptical" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
extern "C" {
  void lidort_fixed_linoptical_c_alloc_init(struct lidort_fixed_linoptical *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_linoptical_c_init_only(struct lidort_fixed_linoptical *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_linoptical_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_linoptical_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_linoptical {
  double* ts_l_deltau_vert_input_;
  int ts_l_deltau_vert_input__f_shapes[2];
  int ts_l_deltau_vert_input__f_byte_size;

  double* ts_l_omega_total_input_;
  int ts_l_omega_total_input__f_shapes[2];
  int ts_l_omega_total_input__f_byte_size;

  double* ts_l_phasfunc_input_dn_;
  int ts_l_phasfunc_input_dn__f_shapes[3];
  int ts_l_phasfunc_input_dn__f_byte_size;

  double* ts_l_phasfunc_input_up_;
  int ts_l_phasfunc_input_up__f_shapes[3];
  int ts_l_phasfunc_input_up__f_byte_size;

  double* ts_l_phasmoms_total_input_;
  int ts_l_phasmoms_total_input__f_shapes[3];
  int ts_l_phasmoms_total_input__f_byte_size;

  
};

// Links to type: "lidort_fixed_linoptical" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
class Lidort_Fixed_Linoptical : public Spurr_Fixed_Linoptical_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Linoptical() : Spurr_Fixed_Linoptical_Base() {
    lidort_fixed_linoptical_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Linoptical(void* allocated_f_type_c) : Spurr_Fixed_Linoptical_Base(allocated_f_type_c) {
    lidort_fixed_linoptical_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Linoptical() {
    if (owns_pointer)
      lidort_fixed_linoptical_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 2>& ts_l_deltau_vert_input() const {
    return ts_l_deltau_vert_input_;
  }

  void ts_l_deltau_vert_input(const blitz::Array<double, 2>& ts_l_deltau_vert_input_in) {
    ts_l_deltau_vert_input_ = ts_l_deltau_vert_input_in;
  }

  
  const blitz::Array<double, 2>& ts_l_omega_total_input() const {
    return ts_l_omega_total_input_;
  }

  void ts_l_omega_total_input(const blitz::Array<double, 2>& ts_l_omega_total_input_in) {
    ts_l_omega_total_input_ = ts_l_omega_total_input_in;
  }

  
  const blitz::Array<double, 3>& ts_l_phasfunc_input_dn() const {
    return ts_l_phasfunc_input_dn_;
  }

  void ts_l_phasfunc_input_dn(const blitz::Array<double, 3>& ts_l_phasfunc_input_dn_in) {
    ts_l_phasfunc_input_dn_ = ts_l_phasfunc_input_dn_in;
  }

  
  const blitz::Array<double, 3>& ts_l_phasfunc_input_up() const {
    return ts_l_phasfunc_input_up_;
  }

  void ts_l_phasfunc_input_up(const blitz::Array<double, 3>& ts_l_phasfunc_input_up_in) {
    ts_l_phasfunc_input_up_ = ts_l_phasfunc_input_up_in;
  }

  
  const blitz::Array<double, 3>& ts_l_phasmoms_total_input() const {
    return ts_l_phasmoms_total_input_;
  }

  void ts_l_phasmoms_total_input(const blitz::Array<double, 3>& ts_l_phasmoms_total_input_in) {
    ts_l_phasmoms_total_input_ = ts_l_phasmoms_total_input_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Linoptical:" << std::endl
      << "   ts_l_deltau_vert_input: " << std::endl << ts_l_deltau_vert_input()  << std::endl
      << "   ts_l_omega_total_input: " << std::endl << ts_l_omega_total_input()  << std::endl
      << "   ts_l_phasfunc_input_dn: " << std::endl << ts_l_phasfunc_input_dn()  << std::endl
      << "   ts_l_phasfunc_input_up: " << std::endl << ts_l_phasfunc_input_up()  << std::endl
      << "ts_l_phasmoms_total_input: " << std::endl << ts_l_phasmoms_total_input()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_l_deltau_vert_input_",sizeof(*transfer_struct_c.ts_l_deltau_vert_input_),transfer_struct_c.ts_l_deltau_vert_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_l_omega_total_input_",sizeof(*transfer_struct_c.ts_l_omega_total_input_),transfer_struct_c.ts_l_omega_total_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_l_phasfunc_input_dn_",sizeof(*transfer_struct_c.ts_l_phasfunc_input_dn_),transfer_struct_c.ts_l_phasfunc_input_dn__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_l_phasfunc_input_up_",sizeof(*transfer_struct_c.ts_l_phasfunc_input_up_),transfer_struct_c.ts_l_phasfunc_input_up__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_l_phasmoms_total_input_",sizeof(*transfer_struct_c.ts_l_phasmoms_total_input_),transfer_struct_c.ts_l_phasmoms_total_input__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_l_deltau_vert_input_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_l_deltau_vert_input_,
      blitz::shape(transfer_struct_c.ts_l_deltau_vert_input__f_shapes[0],
                   transfer_struct_c.ts_l_deltau_vert_input__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_l_omega_total_input_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_l_omega_total_input_,
      blitz::shape(transfer_struct_c.ts_l_omega_total_input__f_shapes[0],
                   transfer_struct_c.ts_l_omega_total_input__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_l_phasfunc_input_dn_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_l_phasfunc_input_dn_,
      blitz::shape(transfer_struct_c.ts_l_phasfunc_input_dn__f_shapes[0],
                   transfer_struct_c.ts_l_phasfunc_input_dn__f_shapes[1],
                   transfer_struct_c.ts_l_phasfunc_input_dn__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_l_phasfunc_input_up_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_l_phasfunc_input_up_,
      blitz::shape(transfer_struct_c.ts_l_phasfunc_input_up__f_shapes[0],
                   transfer_struct_c.ts_l_phasfunc_input_up__f_shapes[1],
                   transfer_struct_c.ts_l_phasfunc_input_up__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_l_phasmoms_total_input_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_l_phasmoms_total_input_,
      blitz::shape(transfer_struct_c.ts_l_phasmoms_total_input__f_shapes[0],
                   transfer_struct_c.ts_l_phasmoms_total_input__f_shapes[1],
                   transfer_struct_c.ts_l_phasmoms_total_input__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_linoptical transfer_struct_c;

  blitz::Array<double, 2> ts_l_deltau_vert_input_;
  blitz::Array<double, 2> ts_l_omega_total_input_;
  blitz::Array<double, 3> ts_l_phasfunc_input_dn_;
  blitz::Array<double, 3> ts_l_phasfunc_input_up_;
  blitz::Array<double, 3> ts_l_phasmoms_total_input_;
  
};

// Links to type: "lidort_fixed_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
extern "C" {
  void lidort_fixed_lininputs_c_alloc_init(struct lidort_fixed_lininputs *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_lininputs_c_init_only(struct lidort_fixed_lininputs *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_lininputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_lininputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_lininputs {
  void* cont_;
  int cont__f_byte_size;

  void* optical_;
  int optical__f_byte_size;

  
};

// Links to type: "lidort_fixed_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
class Lidort_Fixed_Lininputs : public Spurr_Fixed_Lininputs_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Lininputs() : Spurr_Fixed_Lininputs_Base() {
    lidort_fixed_lininputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Lininputs(void* allocated_f_type_c) : Spurr_Fixed_Lininputs_Base(allocated_f_type_c) {
    lidort_fixed_lininputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Lininputs() {
    if (owns_pointer)
      lidort_fixed_lininputs_c_destroy(&fortran_type_c);
  }

  Lidort_Fixed_Lincontrol& cont() {
    return *cont_;
  }

  const Lidort_Fixed_Lincontrol& cont() const {
    return *cont_;
  }

  void cont(Lidort_Fixed_Lincontrol& cont_in) {
    void* src_ptr = cont_in.fortran_type_ptr();
    void* dst_ptr = cont_->fortran_type_ptr();
    lidort_fixed_lincontrol_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Fixed_Linoptical& optical() {
    return *optical_;
  }

  const Lidort_Fixed_Linoptical& optical() const {
    return *optical_;
  }

  void optical(Lidort_Fixed_Linoptical& optical_in) {
    void* src_ptr = optical_in.fortran_type_ptr();
    void* dst_ptr = optical_->fortran_type_ptr();
    lidort_fixed_linoptical_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Lininputs:" << std::endl
      << "   cont: " << cont()  << std::endl
      << "optical: " << optical()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    cont_.reset( new Lidort_Fixed_Lincontrol(transfer_struct_c.cont_) );
    optical_.reset( new Lidort_Fixed_Linoptical(transfer_struct_c.optical_) );
    
  }

  struct lidort_fixed_lininputs transfer_struct_c;

  boost::shared_ptr<Lidort_Fixed_Lincontrol> cont_;
  boost::shared_ptr<Lidort_Fixed_Linoptical> optical_;
  
};

// Links to type: "lidort_modified_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
extern "C" {
  void lidort_modified_lincontrol_c_alloc_init(struct lidort_modified_lincontrol *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_lincontrol_c_init_only(struct lidort_modified_lincontrol *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_lincontrol_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_lincontrol_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_lincontrol {
  int* ts_do_atmos_lbbf_;
  int ts_do_atmos_lbbf__f_byte_size;

  int* ts_do_atmos_linearization_;
  int ts_do_atmos_linearization__f_byte_size;

  int* ts_do_column_linearization_;
  int ts_do_column_linearization__f_byte_size;

  int* ts_do_linearization_;
  int ts_do_linearization__f_byte_size;

  int* ts_do_profile_linearization_;
  int ts_do_profile_linearization__f_byte_size;

  int* ts_do_simulation_only_;
  int ts_do_simulation_only__f_byte_size;

  int* ts_do_sleave_wfs_;
  int ts_do_sleave_wfs__f_byte_size;

  int* ts_do_surface_lbbf_;
  int ts_do_surface_lbbf__f_byte_size;

  int* ts_do_surface_linearization_;
  int ts_do_surface_linearization__f_byte_size;

  
};

// Links to type: "lidort_modified_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
class Lidort_Modified_Lincontrol : public Spurr_Modified_Lincontrol_Base {
public:
  // Allocating constructor
  Lidort_Modified_Lincontrol() : Spurr_Modified_Lincontrol_Base() {
    lidort_modified_lincontrol_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Lincontrol(void* allocated_f_type_c) : Spurr_Modified_Lincontrol_Base(allocated_f_type_c) {
    lidort_modified_lincontrol_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Lincontrol() {
    if (owns_pointer)
      lidort_modified_lincontrol_c_destroy(&fortran_type_c);
  }

  const bool ts_do_atmos_lbbf() const {
    return *transfer_struct_c.ts_do_atmos_lbbf_ != 0;
  }

  void ts_do_atmos_lbbf(const bool& ts_do_atmos_lbbf_in) {
    *transfer_struct_c.ts_do_atmos_lbbf_ = ts_do_atmos_lbbf_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_atmos_linearization() const {
    return *transfer_struct_c.ts_do_atmos_linearization_ != 0;
  }

  void ts_do_atmos_linearization(const bool& ts_do_atmos_linearization_in) {
    *transfer_struct_c.ts_do_atmos_linearization_ = ts_do_atmos_linearization_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_column_linearization() const {
    return *transfer_struct_c.ts_do_column_linearization_ != 0;
  }

  void ts_do_column_linearization(const bool& ts_do_column_linearization_in) {
    *transfer_struct_c.ts_do_column_linearization_ = ts_do_column_linearization_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_linearization() const {
    return *transfer_struct_c.ts_do_linearization_ != 0;
  }

  void ts_do_linearization(const bool& ts_do_linearization_in) {
    *transfer_struct_c.ts_do_linearization_ = ts_do_linearization_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_profile_linearization() const {
    return *transfer_struct_c.ts_do_profile_linearization_ != 0;
  }

  void ts_do_profile_linearization(const bool& ts_do_profile_linearization_in) {
    *transfer_struct_c.ts_do_profile_linearization_ = ts_do_profile_linearization_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_simulation_only() const {
    return *transfer_struct_c.ts_do_simulation_only_ != 0;
  }

  void ts_do_simulation_only(const bool& ts_do_simulation_only_in) {
    *transfer_struct_c.ts_do_simulation_only_ = ts_do_simulation_only_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_sleave_wfs() const {
    return *transfer_struct_c.ts_do_sleave_wfs_ != 0;
  }

  void ts_do_sleave_wfs(const bool& ts_do_sleave_wfs_in) {
    *transfer_struct_c.ts_do_sleave_wfs_ = ts_do_sleave_wfs_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_surface_lbbf() const {
    return *transfer_struct_c.ts_do_surface_lbbf_ != 0;
  }

  void ts_do_surface_lbbf(const bool& ts_do_surface_lbbf_in) {
    *transfer_struct_c.ts_do_surface_lbbf_ = ts_do_surface_lbbf_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_surface_linearization() const {
    return *transfer_struct_c.ts_do_surface_linearization_ != 0;
  }

  void ts_do_surface_linearization(const bool& ts_do_surface_linearization_in) {
    *transfer_struct_c.ts_do_surface_linearization_ = ts_do_surface_linearization_in ? FORTRAN_TRUE_INT : 0;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Lincontrol:" << std::endl
      << "           ts_do_atmos_lbbf: " << ts_do_atmos_lbbf()  << std::endl
      << "  ts_do_atmos_linearization: " << ts_do_atmos_linearization()  << std::endl
      << " ts_do_column_linearization: " << ts_do_column_linearization()  << std::endl
      << "        ts_do_linearization: " << ts_do_linearization()  << std::endl
      << "ts_do_profile_linearization: " << ts_do_profile_linearization()  << std::endl
      << "      ts_do_simulation_only: " << ts_do_simulation_only()  << std::endl
      << "           ts_do_sleave_wfs: " << ts_do_sleave_wfs()  << std::endl
      << "         ts_do_surface_lbbf: " << ts_do_surface_lbbf()  << std::endl
      << "ts_do_surface_linearization: " << ts_do_surface_linearization()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_do_atmos_lbbf_",sizeof(*transfer_struct_c.ts_do_atmos_lbbf_),transfer_struct_c.ts_do_atmos_lbbf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_atmos_linearization_",sizeof(*transfer_struct_c.ts_do_atmos_linearization_),transfer_struct_c.ts_do_atmos_linearization__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_column_linearization_",sizeof(*transfer_struct_c.ts_do_column_linearization_),transfer_struct_c.ts_do_column_linearization__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_linearization_",sizeof(*transfer_struct_c.ts_do_linearization_),transfer_struct_c.ts_do_linearization__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_profile_linearization_",sizeof(*transfer_struct_c.ts_do_profile_linearization_),transfer_struct_c.ts_do_profile_linearization__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_simulation_only_",sizeof(*transfer_struct_c.ts_do_simulation_only_),transfer_struct_c.ts_do_simulation_only__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_sleave_wfs_",sizeof(*transfer_struct_c.ts_do_sleave_wfs_),transfer_struct_c.ts_do_sleave_wfs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_surface_lbbf_",sizeof(*transfer_struct_c.ts_do_surface_lbbf_),transfer_struct_c.ts_do_surface_lbbf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_surface_linearization_",sizeof(*transfer_struct_c.ts_do_surface_linearization_),transfer_struct_c.ts_do_surface_linearization__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_lincontrol transfer_struct_c;

  
};

// Links to type: "lidort_modified_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
extern "C" {
  void lidort_modified_lininputs_c_alloc_init(struct lidort_modified_lininputs *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_lininputs_c_init_only(struct lidort_modified_lininputs *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_lininputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_lininputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_lininputs {
  void* mcont_;
  int mcont__f_byte_size;

  
};

// Links to type: "lidort_modified_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.f90"
class Lidort_Modified_Lininputs : public Spurr_Modified_Lininputs_Base {
public:
  // Allocating constructor
  Lidort_Modified_Lininputs() : Spurr_Modified_Lininputs_Base() {
    lidort_modified_lininputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Lininputs(void* allocated_f_type_c) : Spurr_Modified_Lininputs_Base(allocated_f_type_c) {
    lidort_modified_lininputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Lininputs() {
    if (owns_pointer)
      lidort_modified_lininputs_c_destroy(&fortran_type_c);
  }

  Lidort_Modified_Lincontrol& mcont() {
    return *mcont_;
  }

  const Lidort_Modified_Lincontrol& mcont() const {
    return *mcont_;
  }

  void mcont(Lidort_Modified_Lincontrol& mcont_in) {
    void* src_ptr = mcont_in.fortran_type_ptr();
    void* dst_ptr = mcont_->fortran_type_ptr();
    lidort_modified_lincontrol_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Lininputs:" << std::endl
      << "mcont: " << mcont()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    mcont_.reset( new Lidort_Modified_Lincontrol(transfer_struct_c.mcont_) );
    
  }

  struct lidort_modified_lininputs transfer_struct_c;

  boost::shared_ptr<Lidort_Modified_Lincontrol> mcont_;
  
};

// Links to type: "lidort_linatmos" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.f90"
extern "C" {
  void lidort_linatmos_c_alloc_init(struct lidort_linatmos *transfer_struct_c, void **fortran_type_c);
  void lidort_linatmos_c_init_only(struct lidort_linatmos *transfer_struct_c, void **fortran_type_c);
  void lidort_linatmos_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linatmos_c_destroy(void **fortran_type_c);
  
}

struct lidort_linatmos {
  double* ts_abbwfs_fluxes_;
  int ts_abbwfs_fluxes__f_shapes[4];
  int ts_abbwfs_fluxes__f_byte_size;

  double* ts_abbwfs_jacobians_;
  int ts_abbwfs_jacobians__f_shapes[4];
  int ts_abbwfs_jacobians__f_byte_size;

  double* ts_albmed_fluxes_colwf_;
  int ts_albmed_fluxes_colwf__f_shapes[2];
  int ts_albmed_fluxes_colwf__f_byte_size;

  double* ts_albmed_fluxes_profwf_;
  int ts_albmed_fluxes_profwf__f_shapes[3];
  int ts_albmed_fluxes_profwf__f_byte_size;

  double* ts_albmed_user_colwf_;
  int ts_albmed_user_colwf__f_shapes[2];
  int ts_albmed_user_colwf__f_byte_size;

  double* ts_albmed_user_profwf_;
  int ts_albmed_user_profwf__f_shapes[3];
  int ts_albmed_user_profwf__f_byte_size;

  double* ts_columnwf_;
  int ts_columnwf__f_shapes[4];
  int ts_columnwf__f_byte_size;

  double* ts_dnflux_direct_colwf_;
  int ts_dnflux_direct_colwf__f_shapes[3];
  int ts_dnflux_direct_colwf__f_byte_size;

  double* ts_dnflux_direct_profwf_;
  int ts_dnflux_direct_profwf__f_shapes[4];
  int ts_dnflux_direct_profwf__f_byte_size;

  double* ts_dnmeani_direct_colwf_;
  int ts_dnmeani_direct_colwf__f_shapes[3];
  int ts_dnmeani_direct_colwf__f_byte_size;

  double* ts_dnmeani_direct_profwf_;
  int ts_dnmeani_direct_profwf__f_shapes[4];
  int ts_dnmeani_direct_profwf__f_byte_size;

  double* ts_flux_diffuse_colwf_;
  int ts_flux_diffuse_colwf__f_shapes[4];
  int ts_flux_diffuse_colwf__f_byte_size;

  double* ts_flux_diffuse_profwf_;
  int ts_flux_diffuse_profwf__f_shapes[5];
  int ts_flux_diffuse_profwf__f_byte_size;

  double* ts_lc_layer_mssts_;
  int ts_lc_layer_mssts__f_shapes[3];
  int ts_lc_layer_mssts__f_byte_size;

  double* ts_lc_lostrans_;
  int ts_lc_lostrans__f_shapes[3];
  int ts_lc_lostrans__f_byte_size;

  double* ts_lc_surf_mssts_;
  int ts_lc_surf_mssts__f_shapes[2];
  int ts_lc_surf_mssts__f_byte_size;

  double* ts_lp_layer_mssts_;
  int ts_lp_layer_mssts__f_shapes[4];
  int ts_lp_layer_mssts__f_byte_size;

  double* ts_lp_lostrans_;
  int ts_lp_lostrans__f_shapes[3];
  int ts_lp_lostrans__f_byte_size;

  double* ts_lp_surf_mssts_;
  int ts_lp_surf_mssts__f_shapes[3];
  int ts_lp_surf_mssts__f_byte_size;

  double* ts_meani_diffuse_colwf_;
  int ts_meani_diffuse_colwf__f_shapes[4];
  int ts_meani_diffuse_colwf__f_byte_size;

  double* ts_meani_diffuse_profwf_;
  int ts_meani_diffuse_profwf__f_shapes[5];
  int ts_meani_diffuse_profwf__f_byte_size;

  double* ts_planetary_sbterm_colwf_;
  int ts_planetary_sbterm_colwf__f_shapes[1];
  int ts_planetary_sbterm_colwf__f_byte_size;

  double* ts_planetary_sbterm_profwf_;
  int ts_planetary_sbterm_profwf__f_shapes[2];
  int ts_planetary_sbterm_profwf__f_byte_size;

  double* ts_planetary_transterm_colwf_;
  int ts_planetary_transterm_colwf__f_shapes[2];
  int ts_planetary_transterm_colwf__f_byte_size;

  double* ts_planetary_transterm_profwf_;
  int ts_planetary_transterm_profwf__f_shapes[3];
  int ts_planetary_transterm_profwf__f_byte_size;

  double* ts_profilewf_;
  int ts_profilewf__f_shapes[5];
  int ts_profilewf__f_byte_size;

  double* ts_transbeam_colwf_;
  int ts_transbeam_colwf__f_shapes[2];
  int ts_transbeam_colwf__f_byte_size;

  double* ts_transbeam_profwf_;
  int ts_transbeam_profwf__f_shapes[3];
  int ts_transbeam_profwf__f_byte_size;

  double* ts_trnmed_fluxes_colwf_;
  int ts_trnmed_fluxes_colwf__f_shapes[2];
  int ts_trnmed_fluxes_colwf__f_byte_size;

  double* ts_trnmed_fluxes_profwf_;
  int ts_trnmed_fluxes_profwf__f_shapes[3];
  int ts_trnmed_fluxes_profwf__f_byte_size;

  double* ts_trnmed_user_colwf_;
  int ts_trnmed_user_colwf__f_shapes[2];
  int ts_trnmed_user_colwf__f_byte_size;

  double* ts_trnmed_user_profwf_;
  int ts_trnmed_user_profwf__f_shapes[3];
  int ts_trnmed_user_profwf__f_byte_size;

  
};

// Links to type: "lidort_linatmos" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.f90"
class Lidort_Linatmos : public Spurr_Linatmos_Base {
public:
  // Allocating constructor
  Lidort_Linatmos() : Spurr_Linatmos_Base() {
    lidort_linatmos_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linatmos(void* allocated_f_type_c) : Spurr_Linatmos_Base(allocated_f_type_c) {
    lidort_linatmos_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linatmos() {
    if (owns_pointer)
      lidort_linatmos_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& ts_abbwfs_fluxes() const {
    return ts_abbwfs_fluxes_;
  }

  void ts_abbwfs_fluxes(const blitz::Array<double, 4>& ts_abbwfs_fluxes_in) {
    ts_abbwfs_fluxes_ = ts_abbwfs_fluxes_in;
  }

  
  const blitz::Array<double, 4>& ts_abbwfs_jacobians() const {
    return ts_abbwfs_jacobians_;
  }

  void ts_abbwfs_jacobians(const blitz::Array<double, 4>& ts_abbwfs_jacobians_in) {
    ts_abbwfs_jacobians_ = ts_abbwfs_jacobians_in;
  }

  
  const blitz::Array<double, 2>& ts_albmed_fluxes_colwf() const {
    return ts_albmed_fluxes_colwf_;
  }

  void ts_albmed_fluxes_colwf(const blitz::Array<double, 2>& ts_albmed_fluxes_colwf_in) {
    ts_albmed_fluxes_colwf_ = ts_albmed_fluxes_colwf_in;
  }

  
  const blitz::Array<double, 3>& ts_albmed_fluxes_profwf() const {
    return ts_albmed_fluxes_profwf_;
  }

  void ts_albmed_fluxes_profwf(const blitz::Array<double, 3>& ts_albmed_fluxes_profwf_in) {
    ts_albmed_fluxes_profwf_ = ts_albmed_fluxes_profwf_in;
  }

  
  const blitz::Array<double, 2>& ts_albmed_user_colwf() const {
    return ts_albmed_user_colwf_;
  }

  void ts_albmed_user_colwf(const blitz::Array<double, 2>& ts_albmed_user_colwf_in) {
    ts_albmed_user_colwf_ = ts_albmed_user_colwf_in;
  }

  
  const blitz::Array<double, 3>& ts_albmed_user_profwf() const {
    return ts_albmed_user_profwf_;
  }

  void ts_albmed_user_profwf(const blitz::Array<double, 3>& ts_albmed_user_profwf_in) {
    ts_albmed_user_profwf_ = ts_albmed_user_profwf_in;
  }

  
  const blitz::Array<double, 4>& ts_columnwf() const {
    return ts_columnwf_;
  }

  void ts_columnwf(const blitz::Array<double, 4>& ts_columnwf_in) {
    ts_columnwf_ = ts_columnwf_in;
  }

  
  const blitz::Array<double, 3>& ts_dnflux_direct_colwf() const {
    return ts_dnflux_direct_colwf_;
  }

  void ts_dnflux_direct_colwf(const blitz::Array<double, 3>& ts_dnflux_direct_colwf_in) {
    ts_dnflux_direct_colwf_ = ts_dnflux_direct_colwf_in;
  }

  
  const blitz::Array<double, 4>& ts_dnflux_direct_profwf() const {
    return ts_dnflux_direct_profwf_;
  }

  void ts_dnflux_direct_profwf(const blitz::Array<double, 4>& ts_dnflux_direct_profwf_in) {
    ts_dnflux_direct_profwf_ = ts_dnflux_direct_profwf_in;
  }

  
  const blitz::Array<double, 3>& ts_dnmeani_direct_colwf() const {
    return ts_dnmeani_direct_colwf_;
  }

  void ts_dnmeani_direct_colwf(const blitz::Array<double, 3>& ts_dnmeani_direct_colwf_in) {
    ts_dnmeani_direct_colwf_ = ts_dnmeani_direct_colwf_in;
  }

  
  const blitz::Array<double, 4>& ts_dnmeani_direct_profwf() const {
    return ts_dnmeani_direct_profwf_;
  }

  void ts_dnmeani_direct_profwf(const blitz::Array<double, 4>& ts_dnmeani_direct_profwf_in) {
    ts_dnmeani_direct_profwf_ = ts_dnmeani_direct_profwf_in;
  }

  
  const blitz::Array<double, 4>& ts_flux_diffuse_colwf() const {
    return ts_flux_diffuse_colwf_;
  }

  void ts_flux_diffuse_colwf(const blitz::Array<double, 4>& ts_flux_diffuse_colwf_in) {
    ts_flux_diffuse_colwf_ = ts_flux_diffuse_colwf_in;
  }

  
  const blitz::Array<double, 5>& ts_flux_diffuse_profwf() const {
    return ts_flux_diffuse_profwf_;
  }

  void ts_flux_diffuse_profwf(const blitz::Array<double, 5>& ts_flux_diffuse_profwf_in) {
    ts_flux_diffuse_profwf_ = ts_flux_diffuse_profwf_in;
  }

  
  const blitz::Array<double, 3>& ts_lc_layer_mssts() const {
    return ts_lc_layer_mssts_;
  }

  void ts_lc_layer_mssts(const blitz::Array<double, 3>& ts_lc_layer_mssts_in) {
    ts_lc_layer_mssts_ = ts_lc_layer_mssts_in;
  }

  
  const blitz::Array<double, 3>& ts_lc_lostrans() const {
    return ts_lc_lostrans_;
  }

  void ts_lc_lostrans(const blitz::Array<double, 3>& ts_lc_lostrans_in) {
    ts_lc_lostrans_ = ts_lc_lostrans_in;
  }

  
  const blitz::Array<double, 2>& ts_lc_surf_mssts() const {
    return ts_lc_surf_mssts_;
  }

  void ts_lc_surf_mssts(const blitz::Array<double, 2>& ts_lc_surf_mssts_in) {
    ts_lc_surf_mssts_ = ts_lc_surf_mssts_in;
  }

  
  const blitz::Array<double, 4>& ts_lp_layer_mssts() const {
    return ts_lp_layer_mssts_;
  }

  void ts_lp_layer_mssts(const blitz::Array<double, 4>& ts_lp_layer_mssts_in) {
    ts_lp_layer_mssts_ = ts_lp_layer_mssts_in;
  }

  
  const blitz::Array<double, 3>& ts_lp_lostrans() const {
    return ts_lp_lostrans_;
  }

  void ts_lp_lostrans(const blitz::Array<double, 3>& ts_lp_lostrans_in) {
    ts_lp_lostrans_ = ts_lp_lostrans_in;
  }

  
  const blitz::Array<double, 3>& ts_lp_surf_mssts() const {
    return ts_lp_surf_mssts_;
  }

  void ts_lp_surf_mssts(const blitz::Array<double, 3>& ts_lp_surf_mssts_in) {
    ts_lp_surf_mssts_ = ts_lp_surf_mssts_in;
  }

  
  const blitz::Array<double, 4>& ts_meani_diffuse_colwf() const {
    return ts_meani_diffuse_colwf_;
  }

  void ts_meani_diffuse_colwf(const blitz::Array<double, 4>& ts_meani_diffuse_colwf_in) {
    ts_meani_diffuse_colwf_ = ts_meani_diffuse_colwf_in;
  }

  
  const blitz::Array<double, 5>& ts_meani_diffuse_profwf() const {
    return ts_meani_diffuse_profwf_;
  }

  void ts_meani_diffuse_profwf(const blitz::Array<double, 5>& ts_meani_diffuse_profwf_in) {
    ts_meani_diffuse_profwf_ = ts_meani_diffuse_profwf_in;
  }

  
  const blitz::Array<double, 1>& ts_planetary_sbterm_colwf() const {
    return ts_planetary_sbterm_colwf_;
  }

  void ts_planetary_sbterm_colwf(const blitz::Array<double, 1>& ts_planetary_sbterm_colwf_in) {
    ts_planetary_sbterm_colwf_ = ts_planetary_sbterm_colwf_in;
  }

  
  const blitz::Array<double, 2>& ts_planetary_sbterm_profwf() const {
    return ts_planetary_sbterm_profwf_;
  }

  void ts_planetary_sbterm_profwf(const blitz::Array<double, 2>& ts_planetary_sbterm_profwf_in) {
    ts_planetary_sbterm_profwf_ = ts_planetary_sbterm_profwf_in;
  }

  
  const blitz::Array<double, 2>& ts_planetary_transterm_colwf() const {
    return ts_planetary_transterm_colwf_;
  }

  void ts_planetary_transterm_colwf(const blitz::Array<double, 2>& ts_planetary_transterm_colwf_in) {
    ts_planetary_transterm_colwf_ = ts_planetary_transterm_colwf_in;
  }

  
  const blitz::Array<double, 3>& ts_planetary_transterm_profwf() const {
    return ts_planetary_transterm_profwf_;
  }

  void ts_planetary_transterm_profwf(const blitz::Array<double, 3>& ts_planetary_transterm_profwf_in) {
    ts_planetary_transterm_profwf_ = ts_planetary_transterm_profwf_in;
  }

  
  const blitz::Array<double, 5>& ts_profilewf() const {
    return ts_profilewf_;
  }

  void ts_profilewf(const blitz::Array<double, 5>& ts_profilewf_in) {
    ts_profilewf_ = ts_profilewf_in;
  }

  
  const blitz::Array<double, 2>& ts_transbeam_colwf() const {
    return ts_transbeam_colwf_;
  }

  void ts_transbeam_colwf(const blitz::Array<double, 2>& ts_transbeam_colwf_in) {
    ts_transbeam_colwf_ = ts_transbeam_colwf_in;
  }

  
  const blitz::Array<double, 3>& ts_transbeam_profwf() const {
    return ts_transbeam_profwf_;
  }

  void ts_transbeam_profwf(const blitz::Array<double, 3>& ts_transbeam_profwf_in) {
    ts_transbeam_profwf_ = ts_transbeam_profwf_in;
  }

  
  const blitz::Array<double, 2>& ts_trnmed_fluxes_colwf() const {
    return ts_trnmed_fluxes_colwf_;
  }

  void ts_trnmed_fluxes_colwf(const blitz::Array<double, 2>& ts_trnmed_fluxes_colwf_in) {
    ts_trnmed_fluxes_colwf_ = ts_trnmed_fluxes_colwf_in;
  }

  
  const blitz::Array<double, 3>& ts_trnmed_fluxes_profwf() const {
    return ts_trnmed_fluxes_profwf_;
  }

  void ts_trnmed_fluxes_profwf(const blitz::Array<double, 3>& ts_trnmed_fluxes_profwf_in) {
    ts_trnmed_fluxes_profwf_ = ts_trnmed_fluxes_profwf_in;
  }

  
  const blitz::Array<double, 2>& ts_trnmed_user_colwf() const {
    return ts_trnmed_user_colwf_;
  }

  void ts_trnmed_user_colwf(const blitz::Array<double, 2>& ts_trnmed_user_colwf_in) {
    ts_trnmed_user_colwf_ = ts_trnmed_user_colwf_in;
  }

  
  const blitz::Array<double, 3>& ts_trnmed_user_profwf() const {
    return ts_trnmed_user_profwf_;
  }

  void ts_trnmed_user_profwf(const blitz::Array<double, 3>& ts_trnmed_user_profwf_in) {
    ts_trnmed_user_profwf_ = ts_trnmed_user_profwf_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linatmos:" << std::endl
      << "             ts_abbwfs_fluxes: " << std::endl << ts_abbwfs_fluxes()  << std::endl
      << "          ts_abbwfs_jacobians: " << std::endl << ts_abbwfs_jacobians()  << std::endl
      << "       ts_albmed_fluxes_colwf: " << std::endl << ts_albmed_fluxes_colwf()  << std::endl
      << "      ts_albmed_fluxes_profwf: " << std::endl << ts_albmed_fluxes_profwf()  << std::endl
      << "         ts_albmed_user_colwf: " << std::endl << ts_albmed_user_colwf()  << std::endl
      << "        ts_albmed_user_profwf: " << std::endl << ts_albmed_user_profwf()  << std::endl
      << "                  ts_columnwf: " << std::endl << ts_columnwf()  << std::endl
      << "       ts_dnflux_direct_colwf: " << std::endl << ts_dnflux_direct_colwf()  << std::endl
      << "      ts_dnflux_direct_profwf: " << std::endl << ts_dnflux_direct_profwf()  << std::endl
      << "      ts_dnmeani_direct_colwf: " << std::endl << ts_dnmeani_direct_colwf()  << std::endl
      << "     ts_dnmeani_direct_profwf: " << std::endl << ts_dnmeani_direct_profwf()  << std::endl
      << "        ts_flux_diffuse_colwf: " << std::endl << ts_flux_diffuse_colwf()  << std::endl
      << "       ts_flux_diffuse_profwf: " << std::endl << ts_flux_diffuse_profwf()  << std::endl
      << "            ts_lc_layer_mssts: " << std::endl << ts_lc_layer_mssts()  << std::endl
      << "               ts_lc_lostrans: " << std::endl << ts_lc_lostrans()  << std::endl
      << "             ts_lc_surf_mssts: " << std::endl << ts_lc_surf_mssts()  << std::endl
      << "            ts_lp_layer_mssts: " << std::endl << ts_lp_layer_mssts()  << std::endl
      << "               ts_lp_lostrans: " << std::endl << ts_lp_lostrans()  << std::endl
      << "             ts_lp_surf_mssts: " << std::endl << ts_lp_surf_mssts()  << std::endl
      << "       ts_meani_diffuse_colwf: " << std::endl << ts_meani_diffuse_colwf()  << std::endl
      << "      ts_meani_diffuse_profwf: " << std::endl << ts_meani_diffuse_profwf()  << std::endl
      << "    ts_planetary_sbterm_colwf: " << std::endl << ts_planetary_sbterm_colwf()  << std::endl
      << "   ts_planetary_sbterm_profwf: " << std::endl << ts_planetary_sbterm_profwf()  << std::endl
      << " ts_planetary_transterm_colwf: " << std::endl << ts_planetary_transterm_colwf()  << std::endl
      << "ts_planetary_transterm_profwf: " << std::endl << ts_planetary_transterm_profwf()  << std::endl
      << "                 ts_profilewf: " << std::endl << ts_profilewf()  << std::endl
      << "           ts_transbeam_colwf: " << std::endl << ts_transbeam_colwf()  << std::endl
      << "          ts_transbeam_profwf: " << std::endl << ts_transbeam_profwf()  << std::endl
      << "       ts_trnmed_fluxes_colwf: " << std::endl << ts_trnmed_fluxes_colwf()  << std::endl
      << "      ts_trnmed_fluxes_profwf: " << std::endl << ts_trnmed_fluxes_profwf()  << std::endl
      << "         ts_trnmed_user_colwf: " << std::endl << ts_trnmed_user_colwf()  << std::endl
      << "        ts_trnmed_user_profwf: " << std::endl << ts_trnmed_user_profwf()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_abbwfs_fluxes_",sizeof(*transfer_struct_c.ts_abbwfs_fluxes_),transfer_struct_c.ts_abbwfs_fluxes__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_abbwfs_jacobians_",sizeof(*transfer_struct_c.ts_abbwfs_jacobians_),transfer_struct_c.ts_abbwfs_jacobians__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_albmed_fluxes_colwf_",sizeof(*transfer_struct_c.ts_albmed_fluxes_colwf_),transfer_struct_c.ts_albmed_fluxes_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_albmed_fluxes_profwf_",sizeof(*transfer_struct_c.ts_albmed_fluxes_profwf_),transfer_struct_c.ts_albmed_fluxes_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_albmed_user_colwf_",sizeof(*transfer_struct_c.ts_albmed_user_colwf_),transfer_struct_c.ts_albmed_user_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_albmed_user_profwf_",sizeof(*transfer_struct_c.ts_albmed_user_profwf_),transfer_struct_c.ts_albmed_user_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_columnwf_",sizeof(*transfer_struct_c.ts_columnwf_),transfer_struct_c.ts_columnwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_dnflux_direct_colwf_",sizeof(*transfer_struct_c.ts_dnflux_direct_colwf_),transfer_struct_c.ts_dnflux_direct_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_dnflux_direct_profwf_",sizeof(*transfer_struct_c.ts_dnflux_direct_profwf_),transfer_struct_c.ts_dnflux_direct_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_dnmeani_direct_colwf_",sizeof(*transfer_struct_c.ts_dnmeani_direct_colwf_),transfer_struct_c.ts_dnmeani_direct_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_dnmeani_direct_profwf_",sizeof(*transfer_struct_c.ts_dnmeani_direct_profwf_),transfer_struct_c.ts_dnmeani_direct_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_flux_diffuse_colwf_",sizeof(*transfer_struct_c.ts_flux_diffuse_colwf_),transfer_struct_c.ts_flux_diffuse_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_flux_diffuse_profwf_",sizeof(*transfer_struct_c.ts_flux_diffuse_profwf_),transfer_struct_c.ts_flux_diffuse_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lc_layer_mssts_",sizeof(*transfer_struct_c.ts_lc_layer_mssts_),transfer_struct_c.ts_lc_layer_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lc_lostrans_",sizeof(*transfer_struct_c.ts_lc_lostrans_),transfer_struct_c.ts_lc_lostrans__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lc_surf_mssts_",sizeof(*transfer_struct_c.ts_lc_surf_mssts_),transfer_struct_c.ts_lc_surf_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lp_layer_mssts_",sizeof(*transfer_struct_c.ts_lp_layer_mssts_),transfer_struct_c.ts_lp_layer_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lp_lostrans_",sizeof(*transfer_struct_c.ts_lp_lostrans_),transfer_struct_c.ts_lp_lostrans__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lp_surf_mssts_",sizeof(*transfer_struct_c.ts_lp_surf_mssts_),transfer_struct_c.ts_lp_surf_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_meani_diffuse_colwf_",sizeof(*transfer_struct_c.ts_meani_diffuse_colwf_),transfer_struct_c.ts_meani_diffuse_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_meani_diffuse_profwf_",sizeof(*transfer_struct_c.ts_meani_diffuse_profwf_),transfer_struct_c.ts_meani_diffuse_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_planetary_sbterm_colwf_",sizeof(*transfer_struct_c.ts_planetary_sbterm_colwf_),transfer_struct_c.ts_planetary_sbterm_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_planetary_sbterm_profwf_",sizeof(*transfer_struct_c.ts_planetary_sbterm_profwf_),transfer_struct_c.ts_planetary_sbterm_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_planetary_transterm_colwf_",sizeof(*transfer_struct_c.ts_planetary_transterm_colwf_),transfer_struct_c.ts_planetary_transterm_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_planetary_transterm_profwf_",sizeof(*transfer_struct_c.ts_planetary_transterm_profwf_),transfer_struct_c.ts_planetary_transterm_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_profilewf_",sizeof(*transfer_struct_c.ts_profilewf_),transfer_struct_c.ts_profilewf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_transbeam_colwf_",sizeof(*transfer_struct_c.ts_transbeam_colwf_),transfer_struct_c.ts_transbeam_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_transbeam_profwf_",sizeof(*transfer_struct_c.ts_transbeam_profwf_),transfer_struct_c.ts_transbeam_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_trnmed_fluxes_colwf_",sizeof(*transfer_struct_c.ts_trnmed_fluxes_colwf_),transfer_struct_c.ts_trnmed_fluxes_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_trnmed_fluxes_profwf_",sizeof(*transfer_struct_c.ts_trnmed_fluxes_profwf_),transfer_struct_c.ts_trnmed_fluxes_profwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_trnmed_user_colwf_",sizeof(*transfer_struct_c.ts_trnmed_user_colwf_),transfer_struct_c.ts_trnmed_user_colwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_trnmed_user_profwf_",sizeof(*transfer_struct_c.ts_trnmed_user_profwf_),transfer_struct_c.ts_trnmed_user_profwf__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_abbwfs_fluxes_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_abbwfs_fluxes_,
      blitz::shape(transfer_struct_c.ts_abbwfs_fluxes__f_shapes[0],
                   transfer_struct_c.ts_abbwfs_fluxes__f_shapes[1],
                   transfer_struct_c.ts_abbwfs_fluxes__f_shapes[2],
                   transfer_struct_c.ts_abbwfs_fluxes__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_abbwfs_jacobians_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_abbwfs_jacobians_,
      blitz::shape(transfer_struct_c.ts_abbwfs_jacobians__f_shapes[0],
                   transfer_struct_c.ts_abbwfs_jacobians__f_shapes[1],
                   transfer_struct_c.ts_abbwfs_jacobians__f_shapes[2],
                   transfer_struct_c.ts_abbwfs_jacobians__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_albmed_fluxes_colwf_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_albmed_fluxes_colwf_,
      blitz::shape(transfer_struct_c.ts_albmed_fluxes_colwf__f_shapes[0],
                   transfer_struct_c.ts_albmed_fluxes_colwf__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_albmed_fluxes_profwf_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_albmed_fluxes_profwf_,
      blitz::shape(transfer_struct_c.ts_albmed_fluxes_profwf__f_shapes[0],
                   transfer_struct_c.ts_albmed_fluxes_profwf__f_shapes[1],
                   transfer_struct_c.ts_albmed_fluxes_profwf__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_albmed_user_colwf_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_albmed_user_colwf_,
      blitz::shape(transfer_struct_c.ts_albmed_user_colwf__f_shapes[0],
                   transfer_struct_c.ts_albmed_user_colwf__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_albmed_user_profwf_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_albmed_user_profwf_,
      blitz::shape(transfer_struct_c.ts_albmed_user_profwf__f_shapes[0],
                   transfer_struct_c.ts_albmed_user_profwf__f_shapes[1],
                   transfer_struct_c.ts_albmed_user_profwf__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_columnwf_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_columnwf_,
      blitz::shape(transfer_struct_c.ts_columnwf__f_shapes[0],
                   transfer_struct_c.ts_columnwf__f_shapes[1],
                   transfer_struct_c.ts_columnwf__f_shapes[2],
                   transfer_struct_c.ts_columnwf__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_dnflux_direct_colwf_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_dnflux_direct_colwf_,
      blitz::shape(transfer_struct_c.ts_dnflux_direct_colwf__f_shapes[0],
                   transfer_struct_c.ts_dnflux_direct_colwf__f_shapes[1],
                   transfer_struct_c.ts_dnflux_direct_colwf__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_dnflux_direct_profwf_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_dnflux_direct_profwf_,
      blitz::shape(transfer_struct_c.ts_dnflux_direct_profwf__f_shapes[0],
                   transfer_struct_c.ts_dnflux_direct_profwf__f_shapes[1],
                   transfer_struct_c.ts_dnflux_direct_profwf__f_shapes[2],
                   transfer_struct_c.ts_dnflux_direct_profwf__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_dnmeani_direct_colwf_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_dnmeani_direct_colwf_,
      blitz::shape(transfer_struct_c.ts_dnmeani_direct_colwf__f_shapes[0],
                   transfer_struct_c.ts_dnmeani_direct_colwf__f_shapes[1],
                   transfer_struct_c.ts_dnmeani_direct_colwf__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_dnmeani_direct_profwf_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_dnmeani_direct_profwf_,
      blitz::shape(transfer_struct_c.ts_dnmeani_direct_profwf__f_shapes[0],
                   transfer_struct_c.ts_dnmeani_direct_profwf__f_shapes[1],
                   transfer_struct_c.ts_dnmeani_direct_profwf__f_shapes[2],
                   transfer_struct_c.ts_dnmeani_direct_profwf__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_flux_diffuse_colwf_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_flux_diffuse_colwf_,
      blitz::shape(transfer_struct_c.ts_flux_diffuse_colwf__f_shapes[0],
                   transfer_struct_c.ts_flux_diffuse_colwf__f_shapes[1],
                   transfer_struct_c.ts_flux_diffuse_colwf__f_shapes[2],
                   transfer_struct_c.ts_flux_diffuse_colwf__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_flux_diffuse_profwf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_flux_diffuse_profwf_,
      blitz::shape(transfer_struct_c.ts_flux_diffuse_profwf__f_shapes[0],
                   transfer_struct_c.ts_flux_diffuse_profwf__f_shapes[1],
                   transfer_struct_c.ts_flux_diffuse_profwf__f_shapes[2],
                   transfer_struct_c.ts_flux_diffuse_profwf__f_shapes[3],
                   transfer_struct_c.ts_flux_diffuse_profwf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_lc_layer_mssts_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_lc_layer_mssts_,
      blitz::shape(transfer_struct_c.ts_lc_layer_mssts__f_shapes[0],
                   transfer_struct_c.ts_lc_layer_mssts__f_shapes[1],
                   transfer_struct_c.ts_lc_layer_mssts__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_lc_lostrans_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_lc_lostrans_,
      blitz::shape(transfer_struct_c.ts_lc_lostrans__f_shapes[0],
                   transfer_struct_c.ts_lc_lostrans__f_shapes[1],
                   transfer_struct_c.ts_lc_lostrans__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_lc_surf_mssts_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_lc_surf_mssts_,
      blitz::shape(transfer_struct_c.ts_lc_surf_mssts__f_shapes[0],
                   transfer_struct_c.ts_lc_surf_mssts__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_lp_layer_mssts_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_lp_layer_mssts_,
      blitz::shape(transfer_struct_c.ts_lp_layer_mssts__f_shapes[0],
                   transfer_struct_c.ts_lp_layer_mssts__f_shapes[1],
                   transfer_struct_c.ts_lp_layer_mssts__f_shapes[2],
                   transfer_struct_c.ts_lp_layer_mssts__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_lp_lostrans_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_lp_lostrans_,
      blitz::shape(transfer_struct_c.ts_lp_lostrans__f_shapes[0],
                   transfer_struct_c.ts_lp_lostrans__f_shapes[1],
                   transfer_struct_c.ts_lp_lostrans__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_lp_surf_mssts_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_lp_surf_mssts_,
      blitz::shape(transfer_struct_c.ts_lp_surf_mssts__f_shapes[0],
                   transfer_struct_c.ts_lp_surf_mssts__f_shapes[1],
                   transfer_struct_c.ts_lp_surf_mssts__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_meani_diffuse_colwf_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_meani_diffuse_colwf_,
      blitz::shape(transfer_struct_c.ts_meani_diffuse_colwf__f_shapes[0],
                   transfer_struct_c.ts_meani_diffuse_colwf__f_shapes[1],
                   transfer_struct_c.ts_meani_diffuse_colwf__f_shapes[2],
                   transfer_struct_c.ts_meani_diffuse_colwf__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_meani_diffuse_profwf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_meani_diffuse_profwf_,
      blitz::shape(transfer_struct_c.ts_meani_diffuse_profwf__f_shapes[0],
                   transfer_struct_c.ts_meani_diffuse_profwf__f_shapes[1],
                   transfer_struct_c.ts_meani_diffuse_profwf__f_shapes[2],
                   transfer_struct_c.ts_meani_diffuse_profwf__f_shapes[3],
                   transfer_struct_c.ts_meani_diffuse_profwf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_planetary_sbterm_colwf_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_planetary_sbterm_colwf_,
      blitz::shape(transfer_struct_c.ts_planetary_sbterm_colwf__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_planetary_sbterm_profwf_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_planetary_sbterm_profwf_,
      blitz::shape(transfer_struct_c.ts_planetary_sbterm_profwf__f_shapes[0],
                   transfer_struct_c.ts_planetary_sbterm_profwf__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_planetary_transterm_colwf_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_planetary_transterm_colwf_,
      blitz::shape(transfer_struct_c.ts_planetary_transterm_colwf__f_shapes[0],
                   transfer_struct_c.ts_planetary_transterm_colwf__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_planetary_transterm_profwf_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_planetary_transterm_profwf_,
      blitz::shape(transfer_struct_c.ts_planetary_transterm_profwf__f_shapes[0],
                   transfer_struct_c.ts_planetary_transterm_profwf__f_shapes[1],
                   transfer_struct_c.ts_planetary_transterm_profwf__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_profilewf_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_profilewf_,
      blitz::shape(transfer_struct_c.ts_profilewf__f_shapes[0],
                   transfer_struct_c.ts_profilewf__f_shapes[1],
                   transfer_struct_c.ts_profilewf__f_shapes[2],
                   transfer_struct_c.ts_profilewf__f_shapes[3],
                   transfer_struct_c.ts_profilewf__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    ts_transbeam_colwf_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_transbeam_colwf_,
      blitz::shape(transfer_struct_c.ts_transbeam_colwf__f_shapes[0],
                   transfer_struct_c.ts_transbeam_colwf__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_transbeam_profwf_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_transbeam_profwf_,
      blitz::shape(transfer_struct_c.ts_transbeam_profwf__f_shapes[0],
                   transfer_struct_c.ts_transbeam_profwf__f_shapes[1],
                   transfer_struct_c.ts_transbeam_profwf__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_trnmed_fluxes_colwf_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_trnmed_fluxes_colwf_,
      blitz::shape(transfer_struct_c.ts_trnmed_fluxes_colwf__f_shapes[0],
                   transfer_struct_c.ts_trnmed_fluxes_colwf__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_trnmed_fluxes_profwf_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_trnmed_fluxes_profwf_,
      blitz::shape(transfer_struct_c.ts_trnmed_fluxes_profwf__f_shapes[0],
                   transfer_struct_c.ts_trnmed_fluxes_profwf__f_shapes[1],
                   transfer_struct_c.ts_trnmed_fluxes_profwf__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_trnmed_user_colwf_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_trnmed_user_colwf_,
      blitz::shape(transfer_struct_c.ts_trnmed_user_colwf__f_shapes[0],
                   transfer_struct_c.ts_trnmed_user_colwf__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_trnmed_user_profwf_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_trnmed_user_profwf_,
      blitz::shape(transfer_struct_c.ts_trnmed_user_profwf__f_shapes[0],
                   transfer_struct_c.ts_trnmed_user_profwf__f_shapes[1],
                   transfer_struct_c.ts_trnmed_user_profwf__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linatmos transfer_struct_c;

  blitz::Array<double, 4> ts_abbwfs_fluxes_;
  blitz::Array<double, 4> ts_abbwfs_jacobians_;
  blitz::Array<double, 2> ts_albmed_fluxes_colwf_;
  blitz::Array<double, 3> ts_albmed_fluxes_profwf_;
  blitz::Array<double, 2> ts_albmed_user_colwf_;
  blitz::Array<double, 3> ts_albmed_user_profwf_;
  blitz::Array<double, 4> ts_columnwf_;
  blitz::Array<double, 3> ts_dnflux_direct_colwf_;
  blitz::Array<double, 4> ts_dnflux_direct_profwf_;
  blitz::Array<double, 3> ts_dnmeani_direct_colwf_;
  blitz::Array<double, 4> ts_dnmeani_direct_profwf_;
  blitz::Array<double, 4> ts_flux_diffuse_colwf_;
  blitz::Array<double, 5> ts_flux_diffuse_profwf_;
  blitz::Array<double, 3> ts_lc_layer_mssts_;
  blitz::Array<double, 3> ts_lc_lostrans_;
  blitz::Array<double, 2> ts_lc_surf_mssts_;
  blitz::Array<double, 4> ts_lp_layer_mssts_;
  blitz::Array<double, 3> ts_lp_lostrans_;
  blitz::Array<double, 3> ts_lp_surf_mssts_;
  blitz::Array<double, 4> ts_meani_diffuse_colwf_;
  blitz::Array<double, 5> ts_meani_diffuse_profwf_;
  blitz::Array<double, 1> ts_planetary_sbterm_colwf_;
  blitz::Array<double, 2> ts_planetary_sbterm_profwf_;
  blitz::Array<double, 2> ts_planetary_transterm_colwf_;
  blitz::Array<double, 3> ts_planetary_transterm_profwf_;
  blitz::Array<double, 5> ts_profilewf_;
  blitz::Array<double, 2> ts_transbeam_colwf_;
  blitz::Array<double, 3> ts_transbeam_profwf_;
  blitz::Array<double, 2> ts_trnmed_fluxes_colwf_;
  blitz::Array<double, 3> ts_trnmed_fluxes_profwf_;
  blitz::Array<double, 2> ts_trnmed_user_colwf_;
  blitz::Array<double, 3> ts_trnmed_user_profwf_;
  
};

// Links to type: "lidort_linsurf" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.f90"
extern "C" {
  void lidort_linsurf_c_alloc_init(struct lidort_linsurf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsurf_c_init_only(struct lidort_linsurf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsurf_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsurf_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsurf {
  double* ts_flux_diffuse_surfwf_;
  int ts_flux_diffuse_surfwf__f_shapes[4];
  int ts_flux_diffuse_surfwf__f_byte_size;

  double* ts_ls_layer_mssts_;
  int ts_ls_layer_mssts__f_shapes[3];
  int ts_ls_layer_mssts__f_byte_size;

  double* ts_ls_surf_mssts_;
  int ts_ls_surf_mssts__f_shapes[2];
  int ts_ls_surf_mssts__f_byte_size;

  double* ts_meani_diffuse_surfwf_;
  int ts_meani_diffuse_surfwf__f_shapes[4];
  int ts_meani_diffuse_surfwf__f_byte_size;

  double* ts_sbbwfs_fluxes_;
  int ts_sbbwfs_fluxes__f_shapes[3];
  int ts_sbbwfs_fluxes__f_byte_size;

  double* ts_sbbwfs_jacobians_;
  int ts_sbbwfs_jacobians__f_shapes[3];
  int ts_sbbwfs_jacobians__f_byte_size;

  double* ts_surfacewf_;
  int ts_surfacewf__f_shapes[4];
  int ts_surfacewf__f_byte_size;

  
};

// Links to type: "lidort_linsurf" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.f90"
class Lidort_Linsurf : public Spurr_Linsurf_Base {
public:
  // Allocating constructor
  Lidort_Linsurf() : Spurr_Linsurf_Base() {
    lidort_linsurf_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsurf(void* allocated_f_type_c) : Spurr_Linsurf_Base(allocated_f_type_c) {
    lidort_linsurf_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsurf() {
    if (owns_pointer)
      lidort_linsurf_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& ts_flux_diffuse_surfwf() const {
    return ts_flux_diffuse_surfwf_;
  }

  void ts_flux_diffuse_surfwf(const blitz::Array<double, 4>& ts_flux_diffuse_surfwf_in) {
    ts_flux_diffuse_surfwf_ = ts_flux_diffuse_surfwf_in;
  }

  
  const blitz::Array<double, 3>& ts_ls_layer_mssts() const {
    return ts_ls_layer_mssts_;
  }

  void ts_ls_layer_mssts(const blitz::Array<double, 3>& ts_ls_layer_mssts_in) {
    ts_ls_layer_mssts_ = ts_ls_layer_mssts_in;
  }

  
  const blitz::Array<double, 2>& ts_ls_surf_mssts() const {
    return ts_ls_surf_mssts_;
  }

  void ts_ls_surf_mssts(const blitz::Array<double, 2>& ts_ls_surf_mssts_in) {
    ts_ls_surf_mssts_ = ts_ls_surf_mssts_in;
  }

  
  const blitz::Array<double, 4>& ts_meani_diffuse_surfwf() const {
    return ts_meani_diffuse_surfwf_;
  }

  void ts_meani_diffuse_surfwf(const blitz::Array<double, 4>& ts_meani_diffuse_surfwf_in) {
    ts_meani_diffuse_surfwf_ = ts_meani_diffuse_surfwf_in;
  }

  
  const blitz::Array<double, 3>& ts_sbbwfs_fluxes() const {
    return ts_sbbwfs_fluxes_;
  }

  void ts_sbbwfs_fluxes(const blitz::Array<double, 3>& ts_sbbwfs_fluxes_in) {
    ts_sbbwfs_fluxes_ = ts_sbbwfs_fluxes_in;
  }

  
  const blitz::Array<double, 3>& ts_sbbwfs_jacobians() const {
    return ts_sbbwfs_jacobians_;
  }

  void ts_sbbwfs_jacobians(const blitz::Array<double, 3>& ts_sbbwfs_jacobians_in) {
    ts_sbbwfs_jacobians_ = ts_sbbwfs_jacobians_in;
  }

  
  const blitz::Array<double, 4>& ts_surfacewf() const {
    return ts_surfacewf_;
  }

  void ts_surfacewf(const blitz::Array<double, 4>& ts_surfacewf_in) {
    ts_surfacewf_ = ts_surfacewf_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsurf:" << std::endl
      << " ts_flux_diffuse_surfwf: " << std::endl << ts_flux_diffuse_surfwf()  << std::endl
      << "      ts_ls_layer_mssts: " << std::endl << ts_ls_layer_mssts()  << std::endl
      << "       ts_ls_surf_mssts: " << std::endl << ts_ls_surf_mssts()  << std::endl
      << "ts_meani_diffuse_surfwf: " << std::endl << ts_meani_diffuse_surfwf()  << std::endl
      << "       ts_sbbwfs_fluxes: " << std::endl << ts_sbbwfs_fluxes()  << std::endl
      << "    ts_sbbwfs_jacobians: " << std::endl << ts_sbbwfs_jacobians()  << std::endl
      << "           ts_surfacewf: " << std::endl << ts_surfacewf()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_flux_diffuse_surfwf_",sizeof(*transfer_struct_c.ts_flux_diffuse_surfwf_),transfer_struct_c.ts_flux_diffuse_surfwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_layer_mssts_",sizeof(*transfer_struct_c.ts_ls_layer_mssts_),transfer_struct_c.ts_ls_layer_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_surf_mssts_",sizeof(*transfer_struct_c.ts_ls_surf_mssts_),transfer_struct_c.ts_ls_surf_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_meani_diffuse_surfwf_",sizeof(*transfer_struct_c.ts_meani_diffuse_surfwf_),transfer_struct_c.ts_meani_diffuse_surfwf__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_sbbwfs_fluxes_",sizeof(*transfer_struct_c.ts_sbbwfs_fluxes_),transfer_struct_c.ts_sbbwfs_fluxes__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_sbbwfs_jacobians_",sizeof(*transfer_struct_c.ts_sbbwfs_jacobians_),transfer_struct_c.ts_sbbwfs_jacobians__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_surfacewf_",sizeof(*transfer_struct_c.ts_surfacewf_),transfer_struct_c.ts_surfacewf__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_flux_diffuse_surfwf_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_flux_diffuse_surfwf_,
      blitz::shape(transfer_struct_c.ts_flux_diffuse_surfwf__f_shapes[0],
                   transfer_struct_c.ts_flux_diffuse_surfwf__f_shapes[1],
                   transfer_struct_c.ts_flux_diffuse_surfwf__f_shapes[2],
                   transfer_struct_c.ts_flux_diffuse_surfwf__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_layer_mssts_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_ls_layer_mssts_,
      blitz::shape(transfer_struct_c.ts_ls_layer_mssts__f_shapes[0],
                   transfer_struct_c.ts_ls_layer_mssts__f_shapes[1],
                   transfer_struct_c.ts_ls_layer_mssts__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_ls_surf_mssts_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_ls_surf_mssts_,
      blitz::shape(transfer_struct_c.ts_ls_surf_mssts__f_shapes[0],
                   transfer_struct_c.ts_ls_surf_mssts__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_meani_diffuse_surfwf_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_meani_diffuse_surfwf_,
      blitz::shape(transfer_struct_c.ts_meani_diffuse_surfwf__f_shapes[0],
                   transfer_struct_c.ts_meani_diffuse_surfwf__f_shapes[1],
                   transfer_struct_c.ts_meani_diffuse_surfwf__f_shapes[2],
                   transfer_struct_c.ts_meani_diffuse_surfwf__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_sbbwfs_fluxes_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_sbbwfs_fluxes_,
      blitz::shape(transfer_struct_c.ts_sbbwfs_fluxes__f_shapes[0],
                   transfer_struct_c.ts_sbbwfs_fluxes__f_shapes[1],
                   transfer_struct_c.ts_sbbwfs_fluxes__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_sbbwfs_jacobians_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_sbbwfs_jacobians_,
      blitz::shape(transfer_struct_c.ts_sbbwfs_jacobians__f_shapes[0],
                   transfer_struct_c.ts_sbbwfs_jacobians__f_shapes[1],
                   transfer_struct_c.ts_sbbwfs_jacobians__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_surfacewf_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_surfacewf_,
      blitz::shape(transfer_struct_c.ts_surfacewf__f_shapes[0],
                   transfer_struct_c.ts_surfacewf__f_shapes[1],
                   transfer_struct_c.ts_surfacewf__f_shapes[2],
                   transfer_struct_c.ts_surfacewf__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsurf transfer_struct_c;

  blitz::Array<double, 4> ts_flux_diffuse_surfwf_;
  blitz::Array<double, 3> ts_ls_layer_mssts_;
  blitz::Array<double, 2> ts_ls_surf_mssts_;
  blitz::Array<double, 4> ts_meani_diffuse_surfwf_;
  blitz::Array<double, 3> ts_sbbwfs_fluxes_;
  blitz::Array<double, 3> ts_sbbwfs_jacobians_;
  blitz::Array<double, 4> ts_surfacewf_;
  
};

// Links to type: "lidort_linoutputs" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.f90"
extern "C" {
  void lidort_linoutputs_c_alloc_init(struct lidort_linoutputs *transfer_struct_c, void **fortran_type_c);
  void lidort_linoutputs_c_init_only(struct lidort_linoutputs *transfer_struct_c, void **fortran_type_c);
  void lidort_linoutputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linoutputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_linoutputs {
  void* atmos_;
  int atmos__f_byte_size;

  void* surf_;
  int surf__f_byte_size;

  
};

// Links to type: "lidort_linoutputs" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.f90"
class Lidort_Linoutputs : public Spurr_Linoutputs_Base {
public:
  // Allocating constructor
  Lidort_Linoutputs() : Spurr_Linoutputs_Base() {
    lidort_linoutputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linoutputs(void* allocated_f_type_c) : Spurr_Linoutputs_Base(allocated_f_type_c) {
    lidort_linoutputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linoutputs() {
    if (owns_pointer)
      lidort_linoutputs_c_destroy(&fortran_type_c);
  }

  Lidort_Linatmos& atmos() {
    return *atmos_;
  }

  const Lidort_Linatmos& atmos() const {
    return *atmos_;
  }

  void atmos(Lidort_Linatmos& atmos_in) {
    void* src_ptr = atmos_in.fortran_type_ptr();
    void* dst_ptr = atmos_->fortran_type_ptr();
    lidort_linatmos_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Linsurf& surf() {
    return *surf_;
  }

  const Lidort_Linsurf& surf() const {
    return *surf_;
  }

  void surf(Lidort_Linsurf& surf_in) {
    void* src_ptr = surf_in.fortran_type_ptr();
    void* dst_ptr = surf_->fortran_type_ptr();
    lidort_linsurf_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linoutputs:" << std::endl
      << "atmos: " << atmos()  << std::endl
      << " surf: " << surf()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    atmos_.reset( new Lidort_Linatmos(transfer_struct_c.atmos_) );
    surf_.reset( new Lidort_Linsurf(transfer_struct_c.surf_) );
    
  }

  struct lidort_linoutputs transfer_struct_c;

  boost::shared_ptr<Lidort_Linatmos> atmos_;
  boost::shared_ptr<Lidort_Linsurf> surf_;
  
};

// Links to type: "lidort_linsup_brdf" from module: "lidort_lin_sup_brdf_def_m" in file: "lidort_lin_sup_brdf_def.f90"
extern "C" {
  void lidort_linsup_brdf_c_alloc_init(struct lidort_linsup_brdf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_brdf_c_init_only(struct lidort_linsup_brdf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_brdf_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_brdf_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_brdf {
  double* ts_ls_brdf_f_;
  int ts_ls_brdf_f__f_shapes[4];
  int ts_ls_brdf_f__f_byte_size;

  double* ts_ls_brdf_f_0_;
  int ts_ls_brdf_f_0__f_shapes[4];
  int ts_ls_brdf_f_0__f_byte_size;

  double* ts_ls_emissivity_;
  int ts_ls_emissivity__f_shapes[2];
  int ts_ls_emissivity__f_byte_size;

  double* ts_ls_exactdb_brdfunc_;
  int ts_ls_exactdb_brdfunc__f_shapes[4];
  int ts_ls_exactdb_brdfunc__f_byte_size;

  double* ts_ls_user_brdf_f_;
  int ts_ls_user_brdf_f__f_shapes[4];
  int ts_ls_user_brdf_f__f_byte_size;

  double* ts_ls_user_brdf_f_0_;
  int ts_ls_user_brdf_f_0__f_shapes[4];
  int ts_ls_user_brdf_f_0__f_byte_size;

  double* ts_ls_user_emissivity_;
  int ts_ls_user_emissivity__f_shapes[2];
  int ts_ls_user_emissivity__f_byte_size;

  
};

// Links to type: "lidort_linsup_brdf" from module: "lidort_lin_sup_brdf_def_m" in file: "lidort_lin_sup_brdf_def.f90"
class Lidort_Linsup_Brdf : public Spurr_Linsup_Brdf_Base {
public:
  // Allocating constructor
  Lidort_Linsup_Brdf() : Spurr_Linsup_Brdf_Base() {
    lidort_linsup_brdf_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Brdf(void* allocated_f_type_c) : Spurr_Linsup_Brdf_Base(allocated_f_type_c) {
    lidort_linsup_brdf_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Brdf() {
    if (owns_pointer)
      lidort_linsup_brdf_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& ts_ls_brdf_f() const {
    return ts_ls_brdf_f_;
  }

  void ts_ls_brdf_f(const blitz::Array<double, 4>& ts_ls_brdf_f_in) {
    ts_ls_brdf_f_ = ts_ls_brdf_f_in;
  }

  
  const blitz::Array<double, 4>& ts_ls_brdf_f_0() const {
    return ts_ls_brdf_f_0_;
  }

  void ts_ls_brdf_f_0(const blitz::Array<double, 4>& ts_ls_brdf_f_0_in) {
    ts_ls_brdf_f_0_ = ts_ls_brdf_f_0_in;
  }

  
  const blitz::Array<double, 2>& ts_ls_emissivity() const {
    return ts_ls_emissivity_;
  }

  void ts_ls_emissivity(const blitz::Array<double, 2>& ts_ls_emissivity_in) {
    ts_ls_emissivity_ = ts_ls_emissivity_in;
  }

  
  const blitz::Array<double, 4>& ts_ls_exactdb_brdfunc() const {
    return ts_ls_exactdb_brdfunc_;
  }

  void ts_ls_exactdb_brdfunc(const blitz::Array<double, 4>& ts_ls_exactdb_brdfunc_in) {
    ts_ls_exactdb_brdfunc_ = ts_ls_exactdb_brdfunc_in;
  }

  
  const blitz::Array<double, 4>& ts_ls_user_brdf_f() const {
    return ts_ls_user_brdf_f_;
  }

  void ts_ls_user_brdf_f(const blitz::Array<double, 4>& ts_ls_user_brdf_f_in) {
    ts_ls_user_brdf_f_ = ts_ls_user_brdf_f_in;
  }

  
  const blitz::Array<double, 4>& ts_ls_user_brdf_f_0() const {
    return ts_ls_user_brdf_f_0_;
  }

  void ts_ls_user_brdf_f_0(const blitz::Array<double, 4>& ts_ls_user_brdf_f_0_in) {
    ts_ls_user_brdf_f_0_ = ts_ls_user_brdf_f_0_in;
  }

  
  const blitz::Array<double, 2>& ts_ls_user_emissivity() const {
    return ts_ls_user_emissivity_;
  }

  void ts_ls_user_emissivity(const blitz::Array<double, 2>& ts_ls_user_emissivity_in) {
    ts_ls_user_emissivity_ = ts_ls_user_emissivity_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Brdf:" << std::endl
      << "         ts_ls_brdf_f: " << std::endl << ts_ls_brdf_f()  << std::endl
      << "       ts_ls_brdf_f_0: " << std::endl << ts_ls_brdf_f_0()  << std::endl
      << "     ts_ls_emissivity: " << std::endl << ts_ls_emissivity()  << std::endl
      << "ts_ls_exactdb_brdfunc: " << std::endl << ts_ls_exactdb_brdfunc()  << std::endl
      << "    ts_ls_user_brdf_f: " << std::endl << ts_ls_user_brdf_f()  << std::endl
      << "  ts_ls_user_brdf_f_0: " << std::endl << ts_ls_user_brdf_f_0()  << std::endl
      << "ts_ls_user_emissivity: " << std::endl << ts_ls_user_emissivity()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_ls_brdf_f_",sizeof(*transfer_struct_c.ts_ls_brdf_f_),transfer_struct_c.ts_ls_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_brdf_f_0_",sizeof(*transfer_struct_c.ts_ls_brdf_f_0_),transfer_struct_c.ts_ls_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_emissivity_",sizeof(*transfer_struct_c.ts_ls_emissivity_),transfer_struct_c.ts_ls_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_exactdb_brdfunc_",sizeof(*transfer_struct_c.ts_ls_exactdb_brdfunc_),transfer_struct_c.ts_ls_exactdb_brdfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_user_brdf_f_",sizeof(*transfer_struct_c.ts_ls_user_brdf_f_),transfer_struct_c.ts_ls_user_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_user_brdf_f_0_",sizeof(*transfer_struct_c.ts_ls_user_brdf_f_0_),transfer_struct_c.ts_ls_user_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_ls_user_emissivity_",sizeof(*transfer_struct_c.ts_ls_user_emissivity_),transfer_struct_c.ts_ls_user_emissivity__f_byte_size);
    
  }

  void copy_from_sup(Brdf_Linsup_Outputs& supp_obj) { 
    // This is only safe for if the BRDF and LIDORT
    // sup structures are the exact same structure, but with different
    // names. This MUST be reevaluated on for each version upgrade
    void* sup_ptr = supp_obj.fortran_type_ptr();
    lidort_linsup_brdf_c_copy(&sup_ptr, &fortran_type_c);
  }
      

private:
  void link_blitz_arrays() {
    ts_ls_brdf_f_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_brdf_f_,
      blitz::shape(transfer_struct_c.ts_ls_brdf_f__f_shapes[0],
                   transfer_struct_c.ts_ls_brdf_f__f_shapes[1],
                   transfer_struct_c.ts_ls_brdf_f__f_shapes[2],
                   transfer_struct_c.ts_ls_brdf_f__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_brdf_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_brdf_f_0_,
      blitz::shape(transfer_struct_c.ts_ls_brdf_f_0__f_shapes[0],
                   transfer_struct_c.ts_ls_brdf_f_0__f_shapes[1],
                   transfer_struct_c.ts_ls_brdf_f_0__f_shapes[2],
                   transfer_struct_c.ts_ls_brdf_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_emissivity_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_ls_emissivity_,
      blitz::shape(transfer_struct_c.ts_ls_emissivity__f_shapes[0],
                   transfer_struct_c.ts_ls_emissivity__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_ls_exactdb_brdfunc_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_exactdb_brdfunc_,
      blitz::shape(transfer_struct_c.ts_ls_exactdb_brdfunc__f_shapes[0],
                   transfer_struct_c.ts_ls_exactdb_brdfunc__f_shapes[1],
                   transfer_struct_c.ts_ls_exactdb_brdfunc__f_shapes[2],
                   transfer_struct_c.ts_ls_exactdb_brdfunc__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_user_brdf_f_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_user_brdf_f_,
      blitz::shape(transfer_struct_c.ts_ls_user_brdf_f__f_shapes[0],
                   transfer_struct_c.ts_ls_user_brdf_f__f_shapes[1],
                   transfer_struct_c.ts_ls_user_brdf_f__f_shapes[2],
                   transfer_struct_c.ts_ls_user_brdf_f__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_user_brdf_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_ls_user_brdf_f_0_,
      blitz::shape(transfer_struct_c.ts_ls_user_brdf_f_0__f_shapes[0],
                   transfer_struct_c.ts_ls_user_brdf_f_0__f_shapes[1],
                   transfer_struct_c.ts_ls_user_brdf_f_0__f_shapes[2],
                   transfer_struct_c.ts_ls_user_brdf_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_ls_user_emissivity_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_ls_user_emissivity_,
      blitz::shape(transfer_struct_c.ts_ls_user_emissivity__f_shapes[0],
                   transfer_struct_c.ts_ls_user_emissivity__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsup_brdf transfer_struct_c;

  blitz::Array<double, 4> ts_ls_brdf_f_;
  blitz::Array<double, 4> ts_ls_brdf_f_0_;
  blitz::Array<double, 2> ts_ls_emissivity_;
  blitz::Array<double, 4> ts_ls_exactdb_brdfunc_;
  blitz::Array<double, 4> ts_ls_user_brdf_f_;
  blitz::Array<double, 4> ts_ls_user_brdf_f_0_;
  blitz::Array<double, 2> ts_ls_user_emissivity_;
  
};

// Links to type: "lidort_linsup_sleave" from module: "lidort_lin_sup_sleave_def_m" in file: "lidort_lin_sup_sleave_def.f90"
extern "C" {
  void lidort_linsup_sleave_c_alloc_init(struct lidort_linsup_sleave *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_sleave_c_init_only(struct lidort_linsup_sleave *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_sleave_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_sleave_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_sleave {
  double* ts_lssl_slterm_f_0_;
  int ts_lssl_slterm_f_0__f_shapes[4];
  int ts_lssl_slterm_f_0__f_byte_size;

  double* ts_lssl_slterm_isotropic_;
  int ts_lssl_slterm_isotropic__f_shapes[2];
  int ts_lssl_slterm_isotropic__f_byte_size;

  double* ts_lssl_slterm_userangles_;
  int ts_lssl_slterm_userangles__f_shapes[4];
  int ts_lssl_slterm_userangles__f_byte_size;

  double* ts_lssl_user_slterm_f_0_;
  int ts_lssl_user_slterm_f_0__f_shapes[4];
  int ts_lssl_user_slterm_f_0__f_byte_size;

  
};

// Links to type: "lidort_linsup_sleave" from module: "lidort_lin_sup_sleave_def_m" in file: "lidort_lin_sup_sleave_def.f90"
class Lidort_Linsup_Sleave : public Spurr_Linsup_Sleave_Base {
public:
  // Allocating constructor
  Lidort_Linsup_Sleave() : Spurr_Linsup_Sleave_Base() {
    lidort_linsup_sleave_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Sleave(void* allocated_f_type_c) : Spurr_Linsup_Sleave_Base(allocated_f_type_c) {
    lidort_linsup_sleave_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Sleave() {
    if (owns_pointer)
      lidort_linsup_sleave_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 4>& ts_lssl_slterm_f_0() const {
    return ts_lssl_slterm_f_0_;
  }

  void ts_lssl_slterm_f_0(const blitz::Array<double, 4>& ts_lssl_slterm_f_0_in) {
    ts_lssl_slterm_f_0_ = ts_lssl_slterm_f_0_in;
  }

  
  const blitz::Array<double, 2>& ts_lssl_slterm_isotropic() const {
    return ts_lssl_slterm_isotropic_;
  }

  void ts_lssl_slterm_isotropic(const blitz::Array<double, 2>& ts_lssl_slterm_isotropic_in) {
    ts_lssl_slterm_isotropic_ = ts_lssl_slterm_isotropic_in;
  }

  
  const blitz::Array<double, 4>& ts_lssl_slterm_userangles() const {
    return ts_lssl_slterm_userangles_;
  }

  void ts_lssl_slterm_userangles(const blitz::Array<double, 4>& ts_lssl_slterm_userangles_in) {
    ts_lssl_slterm_userangles_ = ts_lssl_slterm_userangles_in;
  }

  
  const blitz::Array<double, 4>& ts_lssl_user_slterm_f_0() const {
    return ts_lssl_user_slterm_f_0_;
  }

  void ts_lssl_user_slterm_f_0(const blitz::Array<double, 4>& ts_lssl_user_slterm_f_0_in) {
    ts_lssl_user_slterm_f_0_ = ts_lssl_user_slterm_f_0_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Sleave:" << std::endl
      << "       ts_lssl_slterm_f_0: " << std::endl << ts_lssl_slterm_f_0()  << std::endl
      << " ts_lssl_slterm_isotropic: " << std::endl << ts_lssl_slterm_isotropic()  << std::endl
      << "ts_lssl_slterm_userangles: " << std::endl << ts_lssl_slterm_userangles()  << std::endl
      << "  ts_lssl_user_slterm_f_0: " << std::endl << ts_lssl_user_slterm_f_0()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_lssl_slterm_f_0_",sizeof(*transfer_struct_c.ts_lssl_slterm_f_0_),transfer_struct_c.ts_lssl_slterm_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lssl_slterm_isotropic_",sizeof(*transfer_struct_c.ts_lssl_slterm_isotropic_),transfer_struct_c.ts_lssl_slterm_isotropic__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lssl_slterm_userangles_",sizeof(*transfer_struct_c.ts_lssl_slterm_userangles_),transfer_struct_c.ts_lssl_slterm_userangles__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lssl_user_slterm_f_0_",sizeof(*transfer_struct_c.ts_lssl_user_slterm_f_0_),transfer_struct_c.ts_lssl_user_slterm_f_0__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_lssl_slterm_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_lssl_slterm_f_0_,
      blitz::shape(transfer_struct_c.ts_lssl_slterm_f_0__f_shapes[0],
                   transfer_struct_c.ts_lssl_slterm_f_0__f_shapes[1],
                   transfer_struct_c.ts_lssl_slterm_f_0__f_shapes[2],
                   transfer_struct_c.ts_lssl_slterm_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_lssl_slterm_isotropic_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_lssl_slterm_isotropic_,
      blitz::shape(transfer_struct_c.ts_lssl_slterm_isotropic__f_shapes[0],
                   transfer_struct_c.ts_lssl_slterm_isotropic__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_lssl_slterm_userangles_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_lssl_slterm_userangles_,
      blitz::shape(transfer_struct_c.ts_lssl_slterm_userangles__f_shapes[0],
                   transfer_struct_c.ts_lssl_slterm_userangles__f_shapes[1],
                   transfer_struct_c.ts_lssl_slterm_userangles__f_shapes[2],
                   transfer_struct_c.ts_lssl_slterm_userangles__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_lssl_user_slterm_f_0_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_lssl_user_slterm_f_0_,
      blitz::shape(transfer_struct_c.ts_lssl_user_slterm_f_0__f_shapes[0],
                   transfer_struct_c.ts_lssl_user_slterm_f_0__f_shapes[1],
                   transfer_struct_c.ts_lssl_user_slterm_f_0__f_shapes[2],
                   transfer_struct_c.ts_lssl_user_slterm_f_0__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsup_sleave transfer_struct_c;

  blitz::Array<double, 4> ts_lssl_slterm_f_0_;
  blitz::Array<double, 2> ts_lssl_slterm_isotropic_;
  blitz::Array<double, 4> ts_lssl_slterm_userangles_;
  blitz::Array<double, 4> ts_lssl_user_slterm_f_0_;
  
};

// Links to type: "lidort_linsup_ss_atmos" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.f90"
extern "C" {
  void lidort_linsup_ss_atmos_c_alloc_init(struct lidort_linsup_ss_atmos *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_atmos_c_init_only(struct lidort_linsup_ss_atmos *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_atmos_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_ss_atmos_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_ss_atmos {
  double* ts_columnwf_db_;
  int ts_columnwf_db__f_shapes[3];
  int ts_columnwf_db__f_byte_size;

  double* ts_columnwf_ss_;
  int ts_columnwf_ss__f_shapes[4];
  int ts_columnwf_ss__f_byte_size;

  double* ts_profilewf_db_;
  int ts_profilewf_db__f_shapes[4];
  int ts_profilewf_db__f_byte_size;

  double* ts_profilewf_ss_;
  int ts_profilewf_ss__f_shapes[5];
  int ts_profilewf_ss__f_byte_size;

  
};

// Links to type: "lidort_linsup_ss_atmos" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.f90"
class Lidort_Linsup_Ss_Atmos : public Spurr_Type_Structure {
public:
  // Allocating constructor
  Lidort_Linsup_Ss_Atmos() : Spurr_Type_Structure() {
    lidort_linsup_ss_atmos_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Ss_Atmos(void* allocated_f_type_c) : Spurr_Type_Structure(allocated_f_type_c) {
    lidort_linsup_ss_atmos_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Ss_Atmos() {
    if (owns_pointer)
      lidort_linsup_ss_atmos_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_columnwf_db() const {
    return ts_columnwf_db_;
  }

  void ts_columnwf_db(const blitz::Array<double, 3>& ts_columnwf_db_in) {
    ts_columnwf_db_ = ts_columnwf_db_in;
  }

  
  const blitz::Array<double, 4>& ts_columnwf_ss() const {
    return ts_columnwf_ss_;
  }

  void ts_columnwf_ss(const blitz::Array<double, 4>& ts_columnwf_ss_in) {
    ts_columnwf_ss_ = ts_columnwf_ss_in;
  }

  
  const blitz::Array<double, 4>& ts_profilewf_db() const {
    return ts_profilewf_db_;
  }

  void ts_profilewf_db(const blitz::Array<double, 4>& ts_profilewf_db_in) {
    ts_profilewf_db_ = ts_profilewf_db_in;
  }

  
  const blitz::Array<double, 5>& ts_profilewf_ss() const {
    return ts_profilewf_ss_;
  }

  void ts_profilewf_ss(const blitz::Array<double, 5>& ts_profilewf_ss_in) {
    ts_profilewf_ss_ = ts_profilewf_ss_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Ss_Atmos:" << std::endl
      << " ts_columnwf_db: " << std::endl << ts_columnwf_db()  << std::endl
      << " ts_columnwf_ss: " << std::endl << ts_columnwf_ss()  << std::endl
      << "ts_profilewf_db: " << std::endl << ts_profilewf_db()  << std::endl
      << "ts_profilewf_ss: " << std::endl << ts_profilewf_ss()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_columnwf_db_",sizeof(*transfer_struct_c.ts_columnwf_db_),transfer_struct_c.ts_columnwf_db__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_columnwf_ss_",sizeof(*transfer_struct_c.ts_columnwf_ss_),transfer_struct_c.ts_columnwf_ss__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_profilewf_db_",sizeof(*transfer_struct_c.ts_profilewf_db_),transfer_struct_c.ts_profilewf_db__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_profilewf_ss_",sizeof(*transfer_struct_c.ts_profilewf_ss_),transfer_struct_c.ts_profilewf_ss__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_columnwf_db_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_columnwf_db_,
      blitz::shape(transfer_struct_c.ts_columnwf_db__f_shapes[0],
                   transfer_struct_c.ts_columnwf_db__f_shapes[1],
                   transfer_struct_c.ts_columnwf_db__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_columnwf_ss_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_columnwf_ss_,
      blitz::shape(transfer_struct_c.ts_columnwf_ss__f_shapes[0],
                   transfer_struct_c.ts_columnwf_ss__f_shapes[1],
                   transfer_struct_c.ts_columnwf_ss__f_shapes[2],
                   transfer_struct_c.ts_columnwf_ss__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_profilewf_db_.reference(blitz::Array<double, 4>(transfer_struct_c.ts_profilewf_db_,
      blitz::shape(transfer_struct_c.ts_profilewf_db__f_shapes[0],
                   transfer_struct_c.ts_profilewf_db__f_shapes[1],
                   transfer_struct_c.ts_profilewf_db__f_shapes[2],
                   transfer_struct_c.ts_profilewf_db__f_shapes[3]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<4>()));
    ts_profilewf_ss_.reference(blitz::Array<double, 5>(transfer_struct_c.ts_profilewf_ss_,
      blitz::shape(transfer_struct_c.ts_profilewf_ss__f_shapes[0],
                   transfer_struct_c.ts_profilewf_ss__f_shapes[1],
                   transfer_struct_c.ts_profilewf_ss__f_shapes[2],
                   transfer_struct_c.ts_profilewf_ss__f_shapes[3],
                   transfer_struct_c.ts_profilewf_ss__f_shapes[4]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<5>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsup_ss_atmos transfer_struct_c;

  blitz::Array<double, 3> ts_columnwf_db_;
  blitz::Array<double, 4> ts_columnwf_ss_;
  blitz::Array<double, 4> ts_profilewf_db_;
  blitz::Array<double, 5> ts_profilewf_ss_;
  
};

// Links to type: "lidort_linsup_ss_surf" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.f90"
extern "C" {
  void lidort_linsup_ss_surf_c_alloc_init(struct lidort_linsup_ss_surf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_surf_c_init_only(struct lidort_linsup_ss_surf *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_surf_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_ss_surf_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_ss_surf {
  double* ts_surfacewf_db_;
  int ts_surfacewf_db__f_shapes[3];
  int ts_surfacewf_db__f_byte_size;

  
};

// Links to type: "lidort_linsup_ss_surf" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.f90"
class Lidort_Linsup_Ss_Surf : public Spurr_Linsup_Ss_Surf_Base {
public:
  // Allocating constructor
  Lidort_Linsup_Ss_Surf() : Spurr_Linsup_Ss_Surf_Base() {
    lidort_linsup_ss_surf_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Ss_Surf(void* allocated_f_type_c) : Spurr_Linsup_Ss_Surf_Base(allocated_f_type_c) {
    lidort_linsup_ss_surf_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Ss_Surf() {
    if (owns_pointer)
      lidort_linsup_ss_surf_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_surfacewf_db() const {
    return ts_surfacewf_db_;
  }

  void ts_surfacewf_db(const blitz::Array<double, 3>& ts_surfacewf_db_in) {
    ts_surfacewf_db_ = ts_surfacewf_db_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Ss_Surf:" << std::endl
      << "ts_surfacewf_db: " << std::endl << ts_surfacewf_db()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_surfacewf_db_",sizeof(*transfer_struct_c.ts_surfacewf_db_),transfer_struct_c.ts_surfacewf_db__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_surfacewf_db_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_surfacewf_db_,
      blitz::shape(transfer_struct_c.ts_surfacewf_db__f_shapes[0],
                   transfer_struct_c.ts_surfacewf_db__f_shapes[1],
                   transfer_struct_c.ts_surfacewf_db__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_linsup_ss_surf transfer_struct_c;

  blitz::Array<double, 3> ts_surfacewf_db_;
  
};

// Links to type: "lidort_linsup_ss" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.f90"
extern "C" {
  void lidort_linsup_ss_c_alloc_init(struct lidort_linsup_ss *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_c_init_only(struct lidort_linsup_ss *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_ss_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_ss_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_ss {
  void* atmos_;
  int atmos__f_byte_size;

  void* surf_;
  int surf__f_byte_size;

  
};

// Links to type: "lidort_linsup_ss" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.f90"
class Lidort_Linsup_Ss : public Spurr_Linsup_Ss_Base {
public:
  // Allocating constructor
  Lidort_Linsup_Ss() : Spurr_Linsup_Ss_Base() {
    lidort_linsup_ss_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Ss(void* allocated_f_type_c) : Spurr_Linsup_Ss_Base(allocated_f_type_c) {
    lidort_linsup_ss_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Ss() {
    if (owns_pointer)
      lidort_linsup_ss_c_destroy(&fortran_type_c);
  }

  Lidort_Linsup_Ss_Atmos& atmos() {
    return *atmos_;
  }

  const Lidort_Linsup_Ss_Atmos& atmos() const {
    return *atmos_;
  }

  void atmos(Lidort_Linsup_Ss_Atmos& atmos_in) {
    void* src_ptr = atmos_in.fortran_type_ptr();
    void* dst_ptr = atmos_->fortran_type_ptr();
    lidort_linsup_ss_atmos_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Linsup_Ss_Surf& surf() {
    return *surf_;
  }

  const Lidort_Linsup_Ss_Surf& surf() const {
    return *surf_;
  }

  void surf(Lidort_Linsup_Ss_Surf& surf_in) {
    void* src_ptr = surf_in.fortran_type_ptr();
    void* dst_ptr = surf_->fortran_type_ptr();
    lidort_linsup_ss_surf_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Ss:" << std::endl
      << "atmos: " << atmos()  << std::endl
      << " surf: " << surf()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    atmos_.reset( new Lidort_Linsup_Ss_Atmos(transfer_struct_c.atmos_) );
    surf_.reset( new Lidort_Linsup_Ss_Surf(transfer_struct_c.surf_) );
    
  }

  struct lidort_linsup_ss transfer_struct_c;

  boost::shared_ptr<Lidort_Linsup_Ss_Atmos> atmos_;
  boost::shared_ptr<Lidort_Linsup_Ss_Surf> surf_;
  
};

// Links to type: "lidort_linsup_inout" from module: "lidort_lin_sup_inout_def_m" in file: "lidort_lin_sup_def.f90"
extern "C" {
  void lidort_linsup_inout_c_alloc_init(struct lidort_linsup_inout *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_inout_c_init_only(struct lidort_linsup_inout *transfer_struct_c, void **fortran_type_c);
  void lidort_linsup_inout_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_linsup_inout_c_destroy(void **fortran_type_c);
  
}

struct lidort_linsup_inout {
  void* brdf_;
  int brdf__f_byte_size;

  void* sleave_;
  int sleave__f_byte_size;

  void* ss_;
  int ss__f_byte_size;

  
};

// Links to type: "lidort_linsup_inout" from module: "lidort_lin_sup_inout_def_m" in file: "lidort_lin_sup_def.f90"
class Lidort_Linsup_Inout : public Spurr_Linsup_Inout_Base {
public:
  // Allocating constructor
  Lidort_Linsup_Inout() : Spurr_Linsup_Inout_Base() {
    lidort_linsup_inout_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Linsup_Inout(void* allocated_f_type_c) : Spurr_Linsup_Inout_Base(allocated_f_type_c) {
    lidort_linsup_inout_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Linsup_Inout() {
    if (owns_pointer)
      lidort_linsup_inout_c_destroy(&fortran_type_c);
  }

  Lidort_Linsup_Brdf& brdf() {
    return *brdf_;
  }

  const Lidort_Linsup_Brdf& brdf() const {
    return *brdf_;
  }

  void brdf(Lidort_Linsup_Brdf& brdf_in) {
    void* src_ptr = brdf_in.fortran_type_ptr();
    void* dst_ptr = brdf_->fortran_type_ptr();
    lidort_linsup_brdf_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Linsup_Sleave& sleave() {
    return *sleave_;
  }

  const Lidort_Linsup_Sleave& sleave() const {
    return *sleave_;
  }

  void sleave(Lidort_Linsup_Sleave& sleave_in) {
    void* src_ptr = sleave_in.fortran_type_ptr();
    void* dst_ptr = sleave_->fortran_type_ptr();
    lidort_linsup_sleave_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Linsup_Ss& ss() {
    return *ss_;
  }

  const Lidort_Linsup_Ss& ss() const {
    return *ss_;
  }

  void ss(Lidort_Linsup_Ss& ss_in) {
    void* src_ptr = ss_in.fortran_type_ptr();
    void* dst_ptr = ss_->fortran_type_ptr();
    lidort_linsup_ss_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Linsup_Inout:" << std::endl
      << "  brdf: " << brdf()  << std::endl
      << "sleave: " << sleave()  << std::endl
      << "    ss: " << ss()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    brdf_.reset( new Lidort_Linsup_Brdf(transfer_struct_c.brdf_) );
    sleave_.reset( new Lidort_Linsup_Sleave(transfer_struct_c.sleave_) );
    ss_.reset( new Lidort_Linsup_Ss(transfer_struct_c.ss_) );
    
  }

  struct lidort_linsup_inout transfer_struct_c;

  boost::shared_ptr<Lidort_Linsup_Brdf> brdf_;
  boost::shared_ptr<Lidort_Linsup_Sleave> sleave_;
  boost::shared_ptr<Lidort_Linsup_Ss> ss_;
  
};

// Links to type: "lidort_main_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
extern "C" {
  void lidort_main_outputs_c_alloc_init(struct lidort_main_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_main_outputs_c_init_only(struct lidort_main_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_main_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_main_outputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_main_outputs {
  double* ts_albmed_fluxes_;
  int ts_albmed_fluxes__f_shapes[1];
  int ts_albmed_fluxes__f_byte_size;

  double* ts_albmed_user_;
  int ts_albmed_user__f_shapes[1];
  int ts_albmed_user__f_byte_size;

  double* ts_contribs_;
  int ts_contribs__f_shapes[2];
  int ts_contribs__f_byte_size;

  double* ts_dnflux_direct_;
  int ts_dnflux_direct__f_shapes[2];
  int ts_dnflux_direct__f_byte_size;

  double* ts_dnmeani_direct_;
  int ts_dnmeani_direct__f_shapes[2];
  int ts_dnmeani_direct__f_byte_size;

  double* ts_flux_diffuse_;
  int ts_flux_diffuse__f_shapes[3];
  int ts_flux_diffuse__f_byte_size;

  int* ts_fourier_saved_;
  int ts_fourier_saved__f_shapes[1];
  int ts_fourier_saved__f_byte_size;

  double* ts_intensity_;
  int ts_intensity__f_shapes[3];
  int ts_intensity__f_byte_size;

  double* ts_layer_mssts_;
  int ts_layer_mssts__f_shapes[2];
  int ts_layer_mssts__f_byte_size;

  double* ts_lostrans_;
  int ts_lostrans__f_shapes[2];
  int ts_lostrans__f_byte_size;

  double* ts_meani_diffuse_;
  int ts_meani_diffuse__f_shapes[3];
  int ts_meani_diffuse__f_byte_size;

  int* ts_n_geometries_;
  int ts_n_geometries__f_byte_size;

  double* ts_pathgeoms_;
  int ts_pathgeoms__f_shapes[2];
  int ts_pathgeoms__f_byte_size;

  double* ts_planetary_sbterm_;
  int ts_planetary_sbterm__f_byte_size;

  double* ts_planetary_transterm_;
  int ts_planetary_transterm__f_shapes[1];
  int ts_planetary_transterm__f_byte_size;

  double* ts_solarbeam_boatrans_;
  int ts_solarbeam_boatrans__f_shapes[1];
  int ts_solarbeam_boatrans__f_byte_size;

  double* ts_spheralb_;
  int ts_spheralb__f_byte_size;

  double* ts_surf_mssts_;
  int ts_surf_mssts__f_shapes[1];
  int ts_surf_mssts__f_byte_size;

  double* ts_trans1_beam_;
  int ts_trans1_beam__f_shapes[1];
  int ts_trans1_beam__f_byte_size;

  double* ts_trans1_user_;
  int ts_trans1_user__f_shapes[1];
  int ts_trans1_user__f_byte_size;

  double* ts_trnmed_fluxes_;
  int ts_trnmed_fluxes__f_shapes[1];
  int ts_trnmed_fluxes__f_byte_size;

  double* ts_trnmed_user_;
  int ts_trnmed_user__f_shapes[1];
  int ts_trnmed_user__f_byte_size;

  
};

// Links to type: "lidort_main_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
class Lidort_Main_Outputs : public Spurr_Main_Outputs_Base {
public:
  // Allocating constructor
  Lidort_Main_Outputs() : Spurr_Main_Outputs_Base() {
    lidort_main_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Main_Outputs(void* allocated_f_type_c) : Spurr_Main_Outputs_Base(allocated_f_type_c) {
    lidort_main_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Main_Outputs() {
    if (owns_pointer)
      lidort_main_outputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 1>& ts_albmed_fluxes() const {
    return ts_albmed_fluxes_;
  }

  void ts_albmed_fluxes(const blitz::Array<double, 1>& ts_albmed_fluxes_in) {
    ts_albmed_fluxes_ = ts_albmed_fluxes_in;
  }

  
  const blitz::Array<double, 1>& ts_albmed_user() const {
    return ts_albmed_user_;
  }

  void ts_albmed_user(const blitz::Array<double, 1>& ts_albmed_user_in) {
    ts_albmed_user_ = ts_albmed_user_in;
  }

  
  const blitz::Array<double, 2>& ts_contribs() const {
    return ts_contribs_;
  }

  void ts_contribs(const blitz::Array<double, 2>& ts_contribs_in) {
    ts_contribs_ = ts_contribs_in;
  }

  
  const blitz::Array<double, 2>& ts_dnflux_direct() const {
    return ts_dnflux_direct_;
  }

  void ts_dnflux_direct(const blitz::Array<double, 2>& ts_dnflux_direct_in) {
    ts_dnflux_direct_ = ts_dnflux_direct_in;
  }

  
  const blitz::Array<double, 2>& ts_dnmeani_direct() const {
    return ts_dnmeani_direct_;
  }

  void ts_dnmeani_direct(const blitz::Array<double, 2>& ts_dnmeani_direct_in) {
    ts_dnmeani_direct_ = ts_dnmeani_direct_in;
  }

  
  const blitz::Array<double, 3>& ts_flux_diffuse() const {
    return ts_flux_diffuse_;
  }

  void ts_flux_diffuse(const blitz::Array<double, 3>& ts_flux_diffuse_in) {
    ts_flux_diffuse_ = ts_flux_diffuse_in;
  }

  
  const blitz::Array<int, 1>& ts_fourier_saved() const {
    return ts_fourier_saved_;
  }

  void ts_fourier_saved(const blitz::Array<int, 1>& ts_fourier_saved_in) {
    ts_fourier_saved_ = ts_fourier_saved_in;
  }

  
  const blitz::Array<double, 3>& ts_intensity() const {
    return ts_intensity_;
  }

  void ts_intensity(const blitz::Array<double, 3>& ts_intensity_in) {
    ts_intensity_ = ts_intensity_in;
  }

  
  const blitz::Array<double, 2>& ts_layer_mssts() const {
    return ts_layer_mssts_;
  }

  void ts_layer_mssts(const blitz::Array<double, 2>& ts_layer_mssts_in) {
    ts_layer_mssts_ = ts_layer_mssts_in;
  }

  
  const blitz::Array<double, 2>& ts_lostrans() const {
    return ts_lostrans_;
  }

  void ts_lostrans(const blitz::Array<double, 2>& ts_lostrans_in) {
    ts_lostrans_ = ts_lostrans_in;
  }

  
  const blitz::Array<double, 3>& ts_meani_diffuse() const {
    return ts_meani_diffuse_;
  }

  void ts_meani_diffuse(const blitz::Array<double, 3>& ts_meani_diffuse_in) {
    ts_meani_diffuse_ = ts_meani_diffuse_in;
  }

  
  const int& ts_n_geometries() const {
    return *transfer_struct_c.ts_n_geometries_;
  }

  void ts_n_geometries(const int& ts_n_geometries_in) {
    *transfer_struct_c.ts_n_geometries_ = ts_n_geometries_in;
  }

  
  const blitz::Array<double, 2>& ts_pathgeoms() const {
    return ts_pathgeoms_;
  }

  void ts_pathgeoms(const blitz::Array<double, 2>& ts_pathgeoms_in) {
    ts_pathgeoms_ = ts_pathgeoms_in;
  }

  
  const double& ts_planetary_sbterm() const {
    return *transfer_struct_c.ts_planetary_sbterm_;
  }

  void ts_planetary_sbterm(const double& ts_planetary_sbterm_in) {
    *transfer_struct_c.ts_planetary_sbterm_ = ts_planetary_sbterm_in;
  }

  
  const blitz::Array<double, 1>& ts_planetary_transterm() const {
    return ts_planetary_transterm_;
  }

  void ts_planetary_transterm(const blitz::Array<double, 1>& ts_planetary_transterm_in) {
    ts_planetary_transterm_ = ts_planetary_transterm_in;
  }

  
  const blitz::Array<double, 1>& ts_solarbeam_boatrans() const {
    return ts_solarbeam_boatrans_;
  }

  void ts_solarbeam_boatrans(const blitz::Array<double, 1>& ts_solarbeam_boatrans_in) {
    ts_solarbeam_boatrans_ = ts_solarbeam_boatrans_in;
  }

  
  const double& ts_spheralb() const {
    return *transfer_struct_c.ts_spheralb_;
  }

  void ts_spheralb(const double& ts_spheralb_in) {
    *transfer_struct_c.ts_spheralb_ = ts_spheralb_in;
  }

  
  const blitz::Array<double, 1>& ts_surf_mssts() const {
    return ts_surf_mssts_;
  }

  void ts_surf_mssts(const blitz::Array<double, 1>& ts_surf_mssts_in) {
    ts_surf_mssts_ = ts_surf_mssts_in;
  }

  
  const blitz::Array<double, 1>& ts_trans1_beam() const {
    return ts_trans1_beam_;
  }

  void ts_trans1_beam(const blitz::Array<double, 1>& ts_trans1_beam_in) {
    ts_trans1_beam_ = ts_trans1_beam_in;
  }

  
  const blitz::Array<double, 1>& ts_trans1_user() const {
    return ts_trans1_user_;
  }

  void ts_trans1_user(const blitz::Array<double, 1>& ts_trans1_user_in) {
    ts_trans1_user_ = ts_trans1_user_in;
  }

  
  const blitz::Array<double, 1>& ts_trnmed_fluxes() const {
    return ts_trnmed_fluxes_;
  }

  void ts_trnmed_fluxes(const blitz::Array<double, 1>& ts_trnmed_fluxes_in) {
    ts_trnmed_fluxes_ = ts_trnmed_fluxes_in;
  }

  
  const blitz::Array<double, 1>& ts_trnmed_user() const {
    return ts_trnmed_user_;
  }

  void ts_trnmed_user(const blitz::Array<double, 1>& ts_trnmed_user_in) {
    ts_trnmed_user_ = ts_trnmed_user_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Main_Outputs:" << std::endl
      << "      ts_albmed_fluxes: " << std::endl << ts_albmed_fluxes()  << std::endl
      << "        ts_albmed_user: " << std::endl << ts_albmed_user()  << std::endl
      << "           ts_contribs: " << std::endl << ts_contribs()  << std::endl
      << "      ts_dnflux_direct: " << std::endl << ts_dnflux_direct()  << std::endl
      << "     ts_dnmeani_direct: " << std::endl << ts_dnmeani_direct()  << std::endl
      << "       ts_flux_diffuse: " << std::endl << ts_flux_diffuse()  << std::endl
      << "      ts_fourier_saved: " << std::endl << ts_fourier_saved()  << std::endl
      << "          ts_intensity: " << std::endl << ts_intensity()  << std::endl
      << "        ts_layer_mssts: " << std::endl << ts_layer_mssts()  << std::endl
      << "           ts_lostrans: " << std::endl << ts_lostrans()  << std::endl
      << "      ts_meani_diffuse: " << std::endl << ts_meani_diffuse()  << std::endl
      << "       ts_n_geometries: " << ts_n_geometries()  << std::endl
      << "          ts_pathgeoms: " << std::endl << ts_pathgeoms()  << std::endl
      << "   ts_planetary_sbterm: " << ts_planetary_sbterm()  << std::endl
      << "ts_planetary_transterm: " << std::endl << ts_planetary_transterm()  << std::endl
      << " ts_solarbeam_boatrans: " << std::endl << ts_solarbeam_boatrans()  << std::endl
      << "           ts_spheralb: " << ts_spheralb()  << std::endl
      << "         ts_surf_mssts: " << std::endl << ts_surf_mssts()  << std::endl
      << "        ts_trans1_beam: " << std::endl << ts_trans1_beam()  << std::endl
      << "        ts_trans1_user: " << std::endl << ts_trans1_user()  << std::endl
      << "      ts_trnmed_fluxes: " << std::endl << ts_trnmed_fluxes()  << std::endl
      << "        ts_trnmed_user: " << std::endl << ts_trnmed_user()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_albmed_fluxes_",sizeof(*transfer_struct_c.ts_albmed_fluxes_),transfer_struct_c.ts_albmed_fluxes__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_albmed_user_",sizeof(*transfer_struct_c.ts_albmed_user_),transfer_struct_c.ts_albmed_user__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_contribs_",sizeof(*transfer_struct_c.ts_contribs_),transfer_struct_c.ts_contribs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_dnflux_direct_",sizeof(*transfer_struct_c.ts_dnflux_direct_),transfer_struct_c.ts_dnflux_direct__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_dnmeani_direct_",sizeof(*transfer_struct_c.ts_dnmeani_direct_),transfer_struct_c.ts_dnmeani_direct__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_flux_diffuse_",sizeof(*transfer_struct_c.ts_flux_diffuse_),transfer_struct_c.ts_flux_diffuse__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_fourier_saved_",sizeof(*transfer_struct_c.ts_fourier_saved_),transfer_struct_c.ts_fourier_saved__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_intensity_",sizeof(*transfer_struct_c.ts_intensity_),transfer_struct_c.ts_intensity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_layer_mssts_",sizeof(*transfer_struct_c.ts_layer_mssts_),transfer_struct_c.ts_layer_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lostrans_",sizeof(*transfer_struct_c.ts_lostrans_),transfer_struct_c.ts_lostrans__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_meani_diffuse_",sizeof(*transfer_struct_c.ts_meani_diffuse_),transfer_struct_c.ts_meani_diffuse__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_geometries_",sizeof(*transfer_struct_c.ts_n_geometries_),transfer_struct_c.ts_n_geometries__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_pathgeoms_",sizeof(*transfer_struct_c.ts_pathgeoms_),transfer_struct_c.ts_pathgeoms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_planetary_sbterm_",sizeof(*transfer_struct_c.ts_planetary_sbterm_),transfer_struct_c.ts_planetary_sbterm__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_planetary_transterm_",sizeof(*transfer_struct_c.ts_planetary_transterm_),transfer_struct_c.ts_planetary_transterm__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_solarbeam_boatrans_",sizeof(*transfer_struct_c.ts_solarbeam_boatrans_),transfer_struct_c.ts_solarbeam_boatrans__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_spheralb_",sizeof(*transfer_struct_c.ts_spheralb_),transfer_struct_c.ts_spheralb__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_surf_mssts_",sizeof(*transfer_struct_c.ts_surf_mssts_),transfer_struct_c.ts_surf_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_trans1_beam_",sizeof(*transfer_struct_c.ts_trans1_beam_),transfer_struct_c.ts_trans1_beam__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_trans1_user_",sizeof(*transfer_struct_c.ts_trans1_user_),transfer_struct_c.ts_trans1_user__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_trnmed_fluxes_",sizeof(*transfer_struct_c.ts_trnmed_fluxes_),transfer_struct_c.ts_trnmed_fluxes__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_trnmed_user_",sizeof(*transfer_struct_c.ts_trnmed_user_),transfer_struct_c.ts_trnmed_user__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_albmed_fluxes_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_albmed_fluxes_,
      blitz::shape(transfer_struct_c.ts_albmed_fluxes__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_albmed_user_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_albmed_user_,
      blitz::shape(transfer_struct_c.ts_albmed_user__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_contribs_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_contribs_,
      blitz::shape(transfer_struct_c.ts_contribs__f_shapes[0],
                   transfer_struct_c.ts_contribs__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_dnflux_direct_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_dnflux_direct_,
      blitz::shape(transfer_struct_c.ts_dnflux_direct__f_shapes[0],
                   transfer_struct_c.ts_dnflux_direct__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_dnmeani_direct_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_dnmeani_direct_,
      blitz::shape(transfer_struct_c.ts_dnmeani_direct__f_shapes[0],
                   transfer_struct_c.ts_dnmeani_direct__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_flux_diffuse_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_flux_diffuse_,
      blitz::shape(transfer_struct_c.ts_flux_diffuse__f_shapes[0],
                   transfer_struct_c.ts_flux_diffuse__f_shapes[1],
                   transfer_struct_c.ts_flux_diffuse__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_fourier_saved_.reference(blitz::Array<int, 1>(transfer_struct_c.ts_fourier_saved_,
      blitz::shape(transfer_struct_c.ts_fourier_saved__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_intensity_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_intensity_,
      blitz::shape(transfer_struct_c.ts_intensity__f_shapes[0],
                   transfer_struct_c.ts_intensity__f_shapes[1],
                   transfer_struct_c.ts_intensity__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_layer_mssts_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_layer_mssts_,
      blitz::shape(transfer_struct_c.ts_layer_mssts__f_shapes[0],
                   transfer_struct_c.ts_layer_mssts__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_lostrans_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_lostrans_,
      blitz::shape(transfer_struct_c.ts_lostrans__f_shapes[0],
                   transfer_struct_c.ts_lostrans__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_meani_diffuse_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_meani_diffuse_,
      blitz::shape(transfer_struct_c.ts_meani_diffuse__f_shapes[0],
                   transfer_struct_c.ts_meani_diffuse__f_shapes[1],
                   transfer_struct_c.ts_meani_diffuse__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_pathgeoms_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_pathgeoms_,
      blitz::shape(transfer_struct_c.ts_pathgeoms__f_shapes[0],
                   transfer_struct_c.ts_pathgeoms__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_planetary_transterm_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_planetary_transterm_,
      blitz::shape(transfer_struct_c.ts_planetary_transterm__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_solarbeam_boatrans_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_solarbeam_boatrans_,
      blitz::shape(transfer_struct_c.ts_solarbeam_boatrans__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_surf_mssts_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_surf_mssts_,
      blitz::shape(transfer_struct_c.ts_surf_mssts__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_trans1_beam_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_trans1_beam_,
      blitz::shape(transfer_struct_c.ts_trans1_beam__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_trans1_user_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_trans1_user_,
      blitz::shape(transfer_struct_c.ts_trans1_user__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_trnmed_fluxes_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_trnmed_fluxes_,
      blitz::shape(transfer_struct_c.ts_trnmed_fluxes__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_trnmed_user_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_trnmed_user_,
      blitz::shape(transfer_struct_c.ts_trnmed_user__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_main_outputs transfer_struct_c;

  blitz::Array<double, 1> ts_albmed_fluxes_;
  blitz::Array<double, 1> ts_albmed_user_;
  blitz::Array<double, 2> ts_contribs_;
  blitz::Array<double, 2> ts_dnflux_direct_;
  blitz::Array<double, 2> ts_dnmeani_direct_;
  blitz::Array<double, 3> ts_flux_diffuse_;
  blitz::Array<int, 1> ts_fourier_saved_;
  blitz::Array<double, 3> ts_intensity_;
  blitz::Array<double, 2> ts_layer_mssts_;
  blitz::Array<double, 2> ts_lostrans_;
  blitz::Array<double, 3> ts_meani_diffuse_;
  blitz::Array<double, 2> ts_pathgeoms_;
  blitz::Array<double, 1> ts_planetary_transterm_;
  blitz::Array<double, 1> ts_solarbeam_boatrans_;
  blitz::Array<double, 1> ts_surf_mssts_;
  blitz::Array<double, 1> ts_trans1_beam_;
  blitz::Array<double, 1> ts_trans1_user_;
  blitz::Array<double, 1> ts_trnmed_fluxes_;
  blitz::Array<double, 1> ts_trnmed_user_;
  
};

// Links to type: "lidort_wladjusted_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
extern "C" {
  void lidort_wladjusted_outputs_c_alloc_init(struct lidort_wladjusted_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_wladjusted_outputs_c_init_only(struct lidort_wladjusted_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_wladjusted_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_wladjusted_outputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_wladjusted_outputs {
  double* ts_wladjusted_direct_;
  int ts_wladjusted_direct__f_shapes[3];
  int ts_wladjusted_direct__f_byte_size;

  double* ts_wladjusted_f_ords_0_;
  int ts_wladjusted_f_ords_0__f_shapes[3];
  int ts_wladjusted_f_ords_0__f_byte_size;

  double* ts_wladjusted_f_user_0_;
  int ts_wladjusted_f_user_0__f_shapes[3];
  int ts_wladjusted_f_user_0__f_byte_size;

  double* ts_wladjusted_isotropic_;
  int ts_wladjusted_isotropic__f_shapes[1];
  int ts_wladjusted_isotropic__f_byte_size;

  
};

// Links to type: "lidort_wladjusted_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
class Lidort_Wladjusted_Outputs : public Spurr_Wladjusted_Outputs_Base {
public:
  // Allocating constructor
  Lidort_Wladjusted_Outputs() : Spurr_Wladjusted_Outputs_Base() {
    lidort_wladjusted_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Wladjusted_Outputs(void* allocated_f_type_c) : Spurr_Wladjusted_Outputs_Base(allocated_f_type_c) {
    lidort_wladjusted_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Wladjusted_Outputs() {
    if (owns_pointer)
      lidort_wladjusted_outputs_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_wladjusted_direct() const {
    return ts_wladjusted_direct_;
  }

  void ts_wladjusted_direct(const blitz::Array<double, 3>& ts_wladjusted_direct_in) {
    ts_wladjusted_direct_ = ts_wladjusted_direct_in;
  }

  
  const blitz::Array<double, 3>& ts_wladjusted_f_ords_0() const {
    return ts_wladjusted_f_ords_0_;
  }

  void ts_wladjusted_f_ords_0(const blitz::Array<double, 3>& ts_wladjusted_f_ords_0_in) {
    ts_wladjusted_f_ords_0_ = ts_wladjusted_f_ords_0_in;
  }

  
  const blitz::Array<double, 3>& ts_wladjusted_f_user_0() const {
    return ts_wladjusted_f_user_0_;
  }

  void ts_wladjusted_f_user_0(const blitz::Array<double, 3>& ts_wladjusted_f_user_0_in) {
    ts_wladjusted_f_user_0_ = ts_wladjusted_f_user_0_in;
  }

  
  const blitz::Array<double, 1>& ts_wladjusted_isotropic() const {
    return ts_wladjusted_isotropic_;
  }

  void ts_wladjusted_isotropic(const blitz::Array<double, 1>& ts_wladjusted_isotropic_in) {
    ts_wladjusted_isotropic_ = ts_wladjusted_isotropic_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Wladjusted_Outputs:" << std::endl
      << "   ts_wladjusted_direct: " << std::endl << ts_wladjusted_direct()  << std::endl
      << " ts_wladjusted_f_ords_0: " << std::endl << ts_wladjusted_f_ords_0()  << std::endl
      << " ts_wladjusted_f_user_0: " << std::endl << ts_wladjusted_f_user_0()  << std::endl
      << "ts_wladjusted_isotropic: " << std::endl << ts_wladjusted_isotropic()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_wladjusted_direct_",sizeof(*transfer_struct_c.ts_wladjusted_direct_),transfer_struct_c.ts_wladjusted_direct__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_wladjusted_f_ords_0_",sizeof(*transfer_struct_c.ts_wladjusted_f_ords_0_),transfer_struct_c.ts_wladjusted_f_ords_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_wladjusted_f_user_0_",sizeof(*transfer_struct_c.ts_wladjusted_f_user_0_),transfer_struct_c.ts_wladjusted_f_user_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_wladjusted_isotropic_",sizeof(*transfer_struct_c.ts_wladjusted_isotropic_),transfer_struct_c.ts_wladjusted_isotropic__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_wladjusted_direct_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_wladjusted_direct_,
      blitz::shape(transfer_struct_c.ts_wladjusted_direct__f_shapes[0],
                   transfer_struct_c.ts_wladjusted_direct__f_shapes[1],
                   transfer_struct_c.ts_wladjusted_direct__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_wladjusted_f_ords_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_wladjusted_f_ords_0_,
      blitz::shape(transfer_struct_c.ts_wladjusted_f_ords_0__f_shapes[0],
                   transfer_struct_c.ts_wladjusted_f_ords_0__f_shapes[1],
                   transfer_struct_c.ts_wladjusted_f_ords_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_wladjusted_f_user_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_wladjusted_f_user_0_,
      blitz::shape(transfer_struct_c.ts_wladjusted_f_user_0__f_shapes[0],
                   transfer_struct_c.ts_wladjusted_f_user_0__f_shapes[1],
                   transfer_struct_c.ts_wladjusted_f_user_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_wladjusted_isotropic_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_wladjusted_isotropic_,
      blitz::shape(transfer_struct_c.ts_wladjusted_isotropic__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_wladjusted_outputs transfer_struct_c;

  blitz::Array<double, 3> ts_wladjusted_direct_;
  blitz::Array<double, 3> ts_wladjusted_f_ords_0_;
  blitz::Array<double, 3> ts_wladjusted_f_user_0_;
  blitz::Array<double, 1> ts_wladjusted_isotropic_;
  
};

// Links to type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
extern "C" {
  void lidort_exception_handling_c_alloc_init(struct lidort_exception_handling *transfer_struct_c, void **fortran_type_c);
  void lidort_exception_handling_c_init_only(struct lidort_exception_handling *transfer_struct_c, void **fortran_type_c);
  void lidort_exception_handling_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_exception_handling_c_destroy(void **fortran_type_c);
  void exception_handling_ts_actions_get(void **fortran_type_c, const int* ts_actions_in_shape_1, const int* ts_actions_in_len, const char* ts_actions_in);
  void exception_handling_ts_checkmessages_get(void **fortran_type_c, const int* ts_checkmessages_in_shape_1, const int* ts_checkmessages_in_len, const char* ts_checkmessages_in);
  void exception_handling_ts_message_get(void **fortran_type_c, const int* ts_message_in_len, const char* ts_message_in);
  void exception_handling_ts_trace_1_get(void **fortran_type_c, const int* ts_trace_1_in_len, const char* ts_trace_1_in);
  void exception_handling_ts_trace_2_get(void **fortran_type_c, const int* ts_trace_2_in_len, const char* ts_trace_2_in);
  void exception_handling_ts_trace_3_get(void **fortran_type_c, const int* ts_trace_3_in_len, const char* ts_trace_3_in);
  
}

struct lidort_exception_handling {
  
  int ts_actions__f_shapes[1];
  int ts_actions__f_len;

  
  int ts_checkmessages__f_shapes[1];
  int ts_checkmessages__f_len;

  
  int ts_message__f_len;

  int* ts_ncheckmessages_;
  int ts_ncheckmessages__f_byte_size;

  int* ts_status_calculation_;
  int ts_status_calculation__f_byte_size;

  int* ts_status_inputcheck_;
  int ts_status_inputcheck__f_byte_size;

  
  int ts_trace_1__f_len;

  
  int ts_trace_2__f_len;

  
  int ts_trace_3__f_len;

  
};

// Links to type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
class Lidort_Exception_Handling : public Spurr_Exception_Handling_Base {
public:
  // Allocating constructor
  Lidort_Exception_Handling() : Spurr_Exception_Handling_Base() {
    lidort_exception_handling_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Exception_Handling(void* allocated_f_type_c) : Spurr_Exception_Handling_Base(allocated_f_type_c) {
    lidort_exception_handling_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Exception_Handling() {
    if (owns_pointer)
      lidort_exception_handling_c_destroy(&fortran_type_c);
  }

  const std::vector< std::string > ts_actions() const {
    std::vector< std::string > ts_actions_ret;
    blitz::Array<char, 2> ts_actions_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_actions__f_shapes[0], transfer_struct_c.ts_actions__f_len+1, blitz::ColumnMajorArray<2>());
    exception_handling_ts_actions_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_actions__f_shapes[0], &transfer_struct_c.ts_actions__f_len, ts_actions_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_actions_lcl.extent(0); dim_0_idx++)
      ts_actions_ret.push_back( std::string(std::string(ts_actions_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_actions_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_actions_ret;
  }

  
  const std::vector< std::string > ts_checkmessages() const {
    std::vector< std::string > ts_checkmessages_ret;
    blitz::Array<char, 2> ts_checkmessages_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_checkmessages__f_shapes[0], transfer_struct_c.ts_checkmessages__f_len+1, blitz::ColumnMajorArray<2>());
    exception_handling_ts_checkmessages_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_checkmessages__f_shapes[0], &transfer_struct_c.ts_checkmessages__f_len, ts_checkmessages_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_checkmessages_lcl.extent(0); dim_0_idx++)
      ts_checkmessages_ret.push_back( std::string(std::string(ts_checkmessages_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_checkmessages_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_checkmessages_ret;
  }

  
  const std::string ts_message() const {
    std::string ts_message_ret;
    blitz::Array<char, 1> ts_message_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_message__f_len+1, blitz::ColumnMajorArray<1>());
    exception_handling_ts_message_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_message__f_len, ts_message_lcl.dataFirst());
    ts_message_ret = ( std::string(std::string(ts_message_lcl(blitz::Range::all()).begin(), ts_message_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_message_ret;
  }

  
  const int& ts_ncheckmessages() const {
    return *transfer_struct_c.ts_ncheckmessages_;
  }

  void ts_ncheckmessages(const int& ts_ncheckmessages_in) {
    *transfer_struct_c.ts_ncheckmessages_ = ts_ncheckmessages_in;
  }

  
  const int& ts_status_calculation() const {
    return *transfer_struct_c.ts_status_calculation_;
  }

  void ts_status_calculation(const int& ts_status_calculation_in) {
    *transfer_struct_c.ts_status_calculation_ = ts_status_calculation_in;
  }

  
  const int& ts_status_inputcheck() const {
    return *transfer_struct_c.ts_status_inputcheck_;
  }

  void ts_status_inputcheck(const int& ts_status_inputcheck_in) {
    *transfer_struct_c.ts_status_inputcheck_ = ts_status_inputcheck_in;
  }

  
  const std::string ts_trace_1() const {
    std::string ts_trace_1_ret;
    blitz::Array<char, 1> ts_trace_1_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_trace_1__f_len+1, blitz::ColumnMajorArray<1>());
    exception_handling_ts_trace_1_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_trace_1__f_len, ts_trace_1_lcl.dataFirst());
    ts_trace_1_ret = ( std::string(std::string(ts_trace_1_lcl(blitz::Range::all()).begin(), ts_trace_1_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_trace_1_ret;
  }

  
  const std::string ts_trace_2() const {
    std::string ts_trace_2_ret;
    blitz::Array<char, 1> ts_trace_2_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_trace_2__f_len+1, blitz::ColumnMajorArray<1>());
    exception_handling_ts_trace_2_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_trace_2__f_len, ts_trace_2_lcl.dataFirst());
    ts_trace_2_ret = ( std::string(std::string(ts_trace_2_lcl(blitz::Range::all()).begin(), ts_trace_2_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_trace_2_ret;
  }

  
  const std::string ts_trace_3() const {
    std::string ts_trace_3_ret;
    blitz::Array<char, 1> ts_trace_3_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_trace_3__f_len+1, blitz::ColumnMajorArray<1>());
    exception_handling_ts_trace_3_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_trace_3__f_len, ts_trace_3_lcl.dataFirst());
    ts_trace_3_ret = ( std::string(std::string(ts_trace_3_lcl(blitz::Range::all()).begin(), ts_trace_3_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_trace_3_ret;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Exception_Handling:" << std::endl
      << "           ts_actions: " << std::endl;
    std::vector< std::string > ts_actions_lcl = ts_actions();
    for(unsigned int idx = 0; idx < ts_actions_lcl.size(); idx++)
      if ( ts_actions_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_actions_lcl[idx] << "\"" << std::endl;
    output_stream
      << "     ts_checkmessages: " << std::endl;
    std::vector< std::string > ts_checkmessages_lcl = ts_checkmessages();
    for(unsigned int idx = 0; idx < ts_checkmessages_lcl.size(); idx++)
      if ( ts_checkmessages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_checkmessages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "           ts_message: " << "\"" << ts_message() << "\"" << std::endl
      << "    ts_ncheckmessages: " << ts_ncheckmessages()  << std::endl
      << "ts_status_calculation: " << ts_status_calculation()  << std::endl
      << " ts_status_inputcheck: " << ts_status_inputcheck()  << std::endl
      << "           ts_trace_1: " << "\"" << ts_trace_1() << "\"" << std::endl
      << "           ts_trace_2: " << "\"" << ts_trace_2() << "\"" << std::endl
      << "           ts_trace_3: " << "\"" << ts_trace_3() << "\"" << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_ncheckmessages_",sizeof(*transfer_struct_c.ts_ncheckmessages_),transfer_struct_c.ts_ncheckmessages__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_status_calculation_",sizeof(*transfer_struct_c.ts_status_calculation_),transfer_struct_c.ts_status_calculation__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_status_inputcheck_",sizeof(*transfer_struct_c.ts_status_inputcheck_),transfer_struct_c.ts_status_inputcheck__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_exception_handling transfer_struct_c;

  
};

// Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
extern "C" {
  void lidort_input_exception_handling_c_alloc_init(struct lidort_input_exception_handling *transfer_struct_c, void **fortran_type_c);
  void lidort_input_exception_handling_c_init_only(struct lidort_input_exception_handling *transfer_struct_c, void **fortran_type_c);
  void lidort_input_exception_handling_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_input_exception_handling_c_destroy(void **fortran_type_c);
  void input_exception_handling_ts_inputactions_get(void **fortran_type_c, const int* ts_inputactions_in_shape_1, const int* ts_inputactions_in_len, const char* ts_inputactions_in);
  void input_exception_handling_ts_inputmessages_get(void **fortran_type_c, const int* ts_inputmessages_in_shape_1, const int* ts_inputmessages_in_len, const char* ts_inputmessages_in);
  
}

struct lidort_input_exception_handling {
  
  int ts_inputactions__f_shapes[1];
  int ts_inputactions__f_len;

  
  int ts_inputmessages__f_shapes[1];
  int ts_inputmessages__f_len;

  int* ts_ninputmessages_;
  int ts_ninputmessages__f_byte_size;

  int* ts_status_inputread_;
  int ts_status_inputread__f_byte_size;

  
};

// Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
class Lidort_Input_Exception_Handling : public Spurr_Input_Exception_Handling_Base {
public:
  // Allocating constructor
  Lidort_Input_Exception_Handling() : Spurr_Input_Exception_Handling_Base() {
    lidort_input_exception_handling_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Input_Exception_Handling(void* allocated_f_type_c) : Spurr_Input_Exception_Handling_Base(allocated_f_type_c) {
    lidort_input_exception_handling_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Input_Exception_Handling() {
    if (owns_pointer)
      lidort_input_exception_handling_c_destroy(&fortran_type_c);
  }

  const std::vector< std::string > ts_inputactions() const {
    std::vector< std::string > ts_inputactions_ret;
    blitz::Array<char, 2> ts_inputactions_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_inputactions__f_shapes[0], transfer_struct_c.ts_inputactions__f_len+1, blitz::ColumnMajorArray<2>());
    input_exception_handling_ts_inputactions_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_inputactions__f_shapes[0], &transfer_struct_c.ts_inputactions__f_len, ts_inputactions_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_inputactions_lcl.extent(0); dim_0_idx++)
      ts_inputactions_ret.push_back( std::string(std::string(ts_inputactions_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_inputactions_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_inputactions_ret;
  }

  
  const std::vector< std::string > ts_inputmessages() const {
    std::vector< std::string > ts_inputmessages_ret;
    blitz::Array<char, 2> ts_inputmessages_lcl = blitz::Array<char, 2>(transfer_struct_c.ts_inputmessages__f_shapes[0], transfer_struct_c.ts_inputmessages__f_len+1, blitz::ColumnMajorArray<2>());
    input_exception_handling_ts_inputmessages_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_inputmessages__f_shapes[0], &transfer_struct_c.ts_inputmessages__f_len, ts_inputmessages_lcl.dataFirst());
    for(int dim_0_idx = 0; dim_0_idx < ts_inputmessages_lcl.extent(0); dim_0_idx++)
      ts_inputmessages_ret.push_back( std::string(std::string(ts_inputmessages_lcl(dim_0_idx, blitz::Range::all()).begin(), ts_inputmessages_lcl(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return ts_inputmessages_ret;
  }

  
  const int& ts_ninputmessages() const {
    return *transfer_struct_c.ts_ninputmessages_;
  }

  void ts_ninputmessages(const int& ts_ninputmessages_in) {
    *transfer_struct_c.ts_ninputmessages_ = ts_ninputmessages_in;
  }

  
  const int& ts_status_inputread() const {
    return *transfer_struct_c.ts_status_inputread_;
  }

  void ts_status_inputread(const int& ts_status_inputread_in) {
    *transfer_struct_c.ts_status_inputread_ = ts_status_inputread_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Input_Exception_Handling:" << std::endl
      << "    ts_inputactions: " << std::endl;
    std::vector< std::string > ts_inputactions_lcl = ts_inputactions();
    for(unsigned int idx = 0; idx < ts_inputactions_lcl.size(); idx++)
      if ( ts_inputactions_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_inputactions_lcl[idx] << "\"" << std::endl;
    output_stream
      << "   ts_inputmessages: " << std::endl;
    std::vector< std::string > ts_inputmessages_lcl = ts_inputmessages();
    for(unsigned int idx = 0; idx < ts_inputmessages_lcl.size(); idx++)
      if ( ts_inputmessages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << ts_inputmessages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "  ts_ninputmessages: " << ts_ninputmessages()  << std::endl
      << "ts_status_inputread: " << ts_status_inputread()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_ninputmessages_",sizeof(*transfer_struct_c.ts_ninputmessages_),transfer_struct_c.ts_ninputmessages__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_status_inputread_",sizeof(*transfer_struct_c.ts_status_inputread_),transfer_struct_c.ts_status_inputread__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_input_exception_handling transfer_struct_c;

  
};

// Links to type: "lidort_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
extern "C" {
  void lidort_outputs_c_alloc_init(struct lidort_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_outputs_c_init_only(struct lidort_outputs *transfer_struct_c, void **fortran_type_c);
  void lidort_outputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_outputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_outputs {
  void* main_;
  int main__f_byte_size;

  void* status_;
  int status__f_byte_size;

  void* wlout_;
  int wlout__f_byte_size;

  
};

// Links to type: "lidort_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.f90"
class Lidort_Outputs : public Spurr_Outputs_Base {
public:
  // Allocating constructor
  Lidort_Outputs() : Spurr_Outputs_Base() {
    lidort_outputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Outputs(void* allocated_f_type_c) : Spurr_Outputs_Base(allocated_f_type_c) {
    lidort_outputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Outputs() {
    if (owns_pointer)
      lidort_outputs_c_destroy(&fortran_type_c);
  }

  Lidort_Main_Outputs& main() {
    return *main_;
  }

  const Lidort_Main_Outputs& main() const {
    return *main_;
  }

  void main(Lidort_Main_Outputs& main_in) {
    void* src_ptr = main_in.fortran_type_ptr();
    void* dst_ptr = main_->fortran_type_ptr();
    lidort_main_outputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Exception_Handling& status() {
    return *status_;
  }

  const Lidort_Exception_Handling& status() const {
    return *status_;
  }

  void status(Lidort_Exception_Handling& status_in) {
    void* src_ptr = status_in.fortran_type_ptr();
    void* dst_ptr = status_->fortran_type_ptr();
    lidort_exception_handling_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Wladjusted_Outputs& wlout() {
    return *wlout_;
  }

  const Lidort_Wladjusted_Outputs& wlout() const {
    return *wlout_;
  }

  void wlout(Lidort_Wladjusted_Outputs& wlout_in) {
    void* src_ptr = wlout_in.fortran_type_ptr();
    void* dst_ptr = wlout_->fortran_type_ptr();
    lidort_wladjusted_outputs_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Outputs:" << std::endl
      << "  main: " << main()  << std::endl
      << "status: " << status()  << std::endl
      << " wlout: " << wlout()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    main_.reset( new Lidort_Main_Outputs(transfer_struct_c.main_) );
    status_.reset( new Lidort_Exception_Handling(transfer_struct_c.status_) );
    wlout_.reset( new Lidort_Wladjusted_Outputs(transfer_struct_c.wlout_) );
    
  }

  struct lidort_outputs transfer_struct_c;

  boost::shared_ptr<Lidort_Main_Outputs> main_;
  boost::shared_ptr<Lidort_Exception_Handling> status_;
  boost::shared_ptr<Lidort_Wladjusted_Outputs> wlout_;
  
};

// Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def_m" in file: "lidort_sup_brdf_def.f90"
extern "C" {
  void lidort_sup_brdf_c_alloc_init(struct lidort_sup_brdf *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_brdf_c_init_only(struct lidort_sup_brdf *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_brdf_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_sup_brdf_c_destroy(void **fortran_type_c);
  
}

struct lidort_sup_brdf {
  double* ts_brdf_f_;
  int ts_brdf_f__f_shapes[3];
  int ts_brdf_f__f_byte_size;

  double* ts_brdf_f_0_;
  int ts_brdf_f_0__f_shapes[3];
  int ts_brdf_f_0__f_byte_size;

  double* ts_emissivity_;
  int ts_emissivity__f_shapes[1];
  int ts_emissivity__f_byte_size;

  double* ts_exactdb_brdfunc_;
  int ts_exactdb_brdfunc__f_shapes[3];
  int ts_exactdb_brdfunc__f_byte_size;

  double* ts_user_brdf_f_;
  int ts_user_brdf_f__f_shapes[3];
  int ts_user_brdf_f__f_byte_size;

  double* ts_user_brdf_f_0_;
  int ts_user_brdf_f_0__f_shapes[3];
  int ts_user_brdf_f_0__f_byte_size;

  double* ts_user_emissivity_;
  int ts_user_emissivity__f_shapes[1];
  int ts_user_emissivity__f_byte_size;

  
};

// Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def_m" in file: "lidort_sup_brdf_def.f90"
class Lidort_Sup_Brdf : public Spurr_Sup_Brdf_Base {
public:
  // Allocating constructor
  Lidort_Sup_Brdf() : Spurr_Sup_Brdf_Base() {
    lidort_sup_brdf_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Sup_Brdf(void* allocated_f_type_c) : Spurr_Sup_Brdf_Base(allocated_f_type_c) {
    lidort_sup_brdf_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Sup_Brdf() {
    if (owns_pointer)
      lidort_sup_brdf_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_brdf_f() const {
    return ts_brdf_f_;
  }

  void ts_brdf_f(const blitz::Array<double, 3>& ts_brdf_f_in) {
    ts_brdf_f_ = ts_brdf_f_in;
  }

  
  const blitz::Array<double, 3>& ts_brdf_f_0() const {
    return ts_brdf_f_0_;
  }

  void ts_brdf_f_0(const blitz::Array<double, 3>& ts_brdf_f_0_in) {
    ts_brdf_f_0_ = ts_brdf_f_0_in;
  }

  
  const blitz::Array<double, 1>& ts_emissivity() const {
    return ts_emissivity_;
  }

  void ts_emissivity(const blitz::Array<double, 1>& ts_emissivity_in) {
    ts_emissivity_ = ts_emissivity_in;
  }

  
  const blitz::Array<double, 3>& ts_exactdb_brdfunc() const {
    return ts_exactdb_brdfunc_;
  }

  void ts_exactdb_brdfunc(const blitz::Array<double, 3>& ts_exactdb_brdfunc_in) {
    ts_exactdb_brdfunc_ = ts_exactdb_brdfunc_in;
  }

  
  const blitz::Array<double, 3>& ts_user_brdf_f() const {
    return ts_user_brdf_f_;
  }

  void ts_user_brdf_f(const blitz::Array<double, 3>& ts_user_brdf_f_in) {
    ts_user_brdf_f_ = ts_user_brdf_f_in;
  }

  
  const blitz::Array<double, 3>& ts_user_brdf_f_0() const {
    return ts_user_brdf_f_0_;
  }

  void ts_user_brdf_f_0(const blitz::Array<double, 3>& ts_user_brdf_f_0_in) {
    ts_user_brdf_f_0_ = ts_user_brdf_f_0_in;
  }

  
  const blitz::Array<double, 1>& ts_user_emissivity() const {
    return ts_user_emissivity_;
  }

  void ts_user_emissivity(const blitz::Array<double, 1>& ts_user_emissivity_in) {
    ts_user_emissivity_ = ts_user_emissivity_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Sup_Brdf:" << std::endl
      << "         ts_brdf_f: " << std::endl << ts_brdf_f()  << std::endl
      << "       ts_brdf_f_0: " << std::endl << ts_brdf_f_0()  << std::endl
      << "     ts_emissivity: " << std::endl << ts_emissivity()  << std::endl
      << "ts_exactdb_brdfunc: " << std::endl << ts_exactdb_brdfunc()  << std::endl
      << "    ts_user_brdf_f: " << std::endl << ts_user_brdf_f()  << std::endl
      << "  ts_user_brdf_f_0: " << std::endl << ts_user_brdf_f_0()  << std::endl
      << "ts_user_emissivity: " << std::endl << ts_user_emissivity()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_brdf_f_",sizeof(*transfer_struct_c.ts_brdf_f_),transfer_struct_c.ts_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_brdf_f_0_",sizeof(*transfer_struct_c.ts_brdf_f_0_),transfer_struct_c.ts_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_emissivity_",sizeof(*transfer_struct_c.ts_emissivity_),transfer_struct_c.ts_emissivity__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_exactdb_brdfunc_",sizeof(*transfer_struct_c.ts_exactdb_brdfunc_),transfer_struct_c.ts_exactdb_brdfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_brdf_f_",sizeof(*transfer_struct_c.ts_user_brdf_f_),transfer_struct_c.ts_user_brdf_f__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_brdf_f_0_",sizeof(*transfer_struct_c.ts_user_brdf_f_0_),transfer_struct_c.ts_user_brdf_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_emissivity_",sizeof(*transfer_struct_c.ts_user_emissivity_),transfer_struct_c.ts_user_emissivity__f_byte_size);
    
  }

  void copy_from_sup(Brdf_Sup_Outputs& supp_obj) { 
    // This is only safe for if the BRDF and LIDORT
    // sup structures are the exact same structure, but with different
    // names. This MUST be reevaluated on for each version upgrade
    void* sup_ptr = supp_obj.fortran_type_ptr();
    lidort_sup_brdf_c_copy(&sup_ptr, &fortran_type_c);
  }
      

private:
  void link_blitz_arrays() {
    ts_brdf_f_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_brdf_f_,
      blitz::shape(transfer_struct_c.ts_brdf_f__f_shapes[0],
                   transfer_struct_c.ts_brdf_f__f_shapes[1],
                   transfer_struct_c.ts_brdf_f__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_brdf_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_brdf_f_0_,
      blitz::shape(transfer_struct_c.ts_brdf_f_0__f_shapes[0],
                   transfer_struct_c.ts_brdf_f_0__f_shapes[1],
                   transfer_struct_c.ts_brdf_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_emissivity_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_emissivity_,
      blitz::shape(transfer_struct_c.ts_emissivity__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_exactdb_brdfunc_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_exactdb_brdfunc_,
      blitz::shape(transfer_struct_c.ts_exactdb_brdfunc__f_shapes[0],
                   transfer_struct_c.ts_exactdb_brdfunc__f_shapes[1],
                   transfer_struct_c.ts_exactdb_brdfunc__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_user_brdf_f_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_user_brdf_f_,
      blitz::shape(transfer_struct_c.ts_user_brdf_f__f_shapes[0],
                   transfer_struct_c.ts_user_brdf_f__f_shapes[1],
                   transfer_struct_c.ts_user_brdf_f__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_user_brdf_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_user_brdf_f_0_,
      blitz::shape(transfer_struct_c.ts_user_brdf_f_0__f_shapes[0],
                   transfer_struct_c.ts_user_brdf_f_0__f_shapes[1],
                   transfer_struct_c.ts_user_brdf_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_user_emissivity_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_user_emissivity_,
      blitz::shape(transfer_struct_c.ts_user_emissivity__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_sup_brdf transfer_struct_c;

  blitz::Array<double, 3> ts_brdf_f_;
  blitz::Array<double, 3> ts_brdf_f_0_;
  blitz::Array<double, 1> ts_emissivity_;
  blitz::Array<double, 3> ts_exactdb_brdfunc_;
  blitz::Array<double, 3> ts_user_brdf_f_;
  blitz::Array<double, 3> ts_user_brdf_f_0_;
  blitz::Array<double, 1> ts_user_emissivity_;
  
};

// Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def_m" in file: "lidort_sup_sleave_def.f90"
extern "C" {
  void lidort_sup_sleave_c_alloc_init(struct lidort_sup_sleave *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_sleave_c_init_only(struct lidort_sup_sleave *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_sleave_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_sup_sleave_c_destroy(void **fortran_type_c);
  
}

struct lidort_sup_sleave {
  double* ts_slterm_f_0_;
  int ts_slterm_f_0__f_shapes[3];
  int ts_slterm_f_0__f_byte_size;

  double* ts_slterm_isotropic_;
  int ts_slterm_isotropic__f_shapes[1];
  int ts_slterm_isotropic__f_byte_size;

  double* ts_slterm_userangles_;
  int ts_slterm_userangles__f_shapes[3];
  int ts_slterm_userangles__f_byte_size;

  double* ts_user_slterm_f_0_;
  int ts_user_slterm_f_0__f_shapes[3];
  int ts_user_slterm_f_0__f_byte_size;

  
};

// Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def_m" in file: "lidort_sup_sleave_def.f90"
class Lidort_Sup_Sleave : public Spurr_Sup_Sleave_Base {
public:
  // Allocating constructor
  Lidort_Sup_Sleave() : Spurr_Sup_Sleave_Base() {
    lidort_sup_sleave_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Sup_Sleave(void* allocated_f_type_c) : Spurr_Sup_Sleave_Base(allocated_f_type_c) {
    lidort_sup_sleave_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Sup_Sleave() {
    if (owns_pointer)
      lidort_sup_sleave_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 3>& ts_slterm_f_0() const {
    return ts_slterm_f_0_;
  }

  void ts_slterm_f_0(const blitz::Array<double, 3>& ts_slterm_f_0_in) {
    ts_slterm_f_0_ = ts_slterm_f_0_in;
  }

  
  const blitz::Array<double, 1>& ts_slterm_isotropic() const {
    return ts_slterm_isotropic_;
  }

  void ts_slterm_isotropic(const blitz::Array<double, 1>& ts_slterm_isotropic_in) {
    ts_slterm_isotropic_ = ts_slterm_isotropic_in;
  }

  
  const blitz::Array<double, 3>& ts_slterm_userangles() const {
    return ts_slterm_userangles_;
  }

  void ts_slterm_userangles(const blitz::Array<double, 3>& ts_slterm_userangles_in) {
    ts_slterm_userangles_ = ts_slterm_userangles_in;
  }

  
  const blitz::Array<double, 3>& ts_user_slterm_f_0() const {
    return ts_user_slterm_f_0_;
  }

  void ts_user_slterm_f_0(const blitz::Array<double, 3>& ts_user_slterm_f_0_in) {
    ts_user_slterm_f_0_ = ts_user_slterm_f_0_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Sup_Sleave:" << std::endl
      << "       ts_slterm_f_0: " << std::endl << ts_slterm_f_0()  << std::endl
      << " ts_slterm_isotropic: " << std::endl << ts_slterm_isotropic()  << std::endl
      << "ts_slterm_userangles: " << std::endl << ts_slterm_userangles()  << std::endl
      << "  ts_user_slterm_f_0: " << std::endl << ts_user_slterm_f_0()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_slterm_f_0_",sizeof(*transfer_struct_c.ts_slterm_f_0_),transfer_struct_c.ts_slterm_f_0__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_slterm_isotropic_",sizeof(*transfer_struct_c.ts_slterm_isotropic_),transfer_struct_c.ts_slterm_isotropic__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_slterm_userangles_",sizeof(*transfer_struct_c.ts_slterm_userangles_),transfer_struct_c.ts_slterm_userangles__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_slterm_f_0_",sizeof(*transfer_struct_c.ts_user_slterm_f_0_),transfer_struct_c.ts_user_slterm_f_0__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_slterm_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_slterm_f_0_,
      blitz::shape(transfer_struct_c.ts_slterm_f_0__f_shapes[0],
                   transfer_struct_c.ts_slterm_f_0__f_shapes[1],
                   transfer_struct_c.ts_slterm_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_slterm_isotropic_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_slterm_isotropic_,
      blitz::shape(transfer_struct_c.ts_slterm_isotropic__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_slterm_userangles_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_slterm_userangles_,
      blitz::shape(transfer_struct_c.ts_slterm_userangles__f_shapes[0],
                   transfer_struct_c.ts_slterm_userangles__f_shapes[1],
                   transfer_struct_c.ts_slterm_userangles__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    ts_user_slterm_f_0_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_user_slterm_f_0_,
      blitz::shape(transfer_struct_c.ts_user_slterm_f_0__f_shapes[0],
                   transfer_struct_c.ts_user_slterm_f_0__f_shapes[1],
                   transfer_struct_c.ts_user_slterm_f_0__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_sup_sleave transfer_struct_c;

  blitz::Array<double, 3> ts_slterm_f_0_;
  blitz::Array<double, 1> ts_slterm_isotropic_;
  blitz::Array<double, 3> ts_slterm_userangles_;
  blitz::Array<double, 3> ts_user_slterm_f_0_;
  
};

// Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def_m" in file: "lidort_sup_ss_def.f90"
extern "C" {
  void lidort_sup_ss_c_alloc_init(struct lidort_sup_ss *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_ss_c_init_only(struct lidort_sup_ss *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_ss_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_sup_ss_c_destroy(void **fortran_type_c);
  
}

struct lidort_sup_ss {
  double* ts_contribs_ss_;
  int ts_contribs_ss__f_shapes[2];
  int ts_contribs_ss__f_byte_size;

  double* ts_intensity_db_;
  int ts_intensity_db__f_shapes[2];
  int ts_intensity_db__f_byte_size;

  double* ts_intensity_ss_;
  int ts_intensity_ss__f_shapes[3];
  int ts_intensity_ss__f_byte_size;

  
};

// Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def_m" in file: "lidort_sup_ss_def.f90"
class Lidort_Sup_Ss : public Spurr_Sup_Ss_Base {
public:
  // Allocating constructor
  Lidort_Sup_Ss() : Spurr_Sup_Ss_Base() {
    lidort_sup_ss_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Sup_Ss(void* allocated_f_type_c) : Spurr_Sup_Ss_Base(allocated_f_type_c) {
    lidort_sup_ss_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Sup_Ss() {
    if (owns_pointer)
      lidort_sup_ss_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 2>& ts_contribs_ss() const {
    return ts_contribs_ss_;
  }

  void ts_contribs_ss(const blitz::Array<double, 2>& ts_contribs_ss_in) {
    ts_contribs_ss_ = ts_contribs_ss_in;
  }

  
  const blitz::Array<double, 2>& ts_intensity_db() const {
    return ts_intensity_db_;
  }

  void ts_intensity_db(const blitz::Array<double, 2>& ts_intensity_db_in) {
    ts_intensity_db_ = ts_intensity_db_in;
  }

  
  const blitz::Array<double, 3>& ts_intensity_ss() const {
    return ts_intensity_ss_;
  }

  void ts_intensity_ss(const blitz::Array<double, 3>& ts_intensity_ss_in) {
    ts_intensity_ss_ = ts_intensity_ss_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Sup_Ss:" << std::endl
      << " ts_contribs_ss: " << std::endl << ts_contribs_ss()  << std::endl
      << "ts_intensity_db: " << std::endl << ts_intensity_db()  << std::endl
      << "ts_intensity_ss: " << std::endl << ts_intensity_ss()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_contribs_ss_",sizeof(*transfer_struct_c.ts_contribs_ss_),transfer_struct_c.ts_contribs_ss__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_intensity_db_",sizeof(*transfer_struct_c.ts_intensity_db_),transfer_struct_c.ts_intensity_db__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_intensity_ss_",sizeof(*transfer_struct_c.ts_intensity_ss_),transfer_struct_c.ts_intensity_ss__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_contribs_ss_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_contribs_ss_,
      blitz::shape(transfer_struct_c.ts_contribs_ss__f_shapes[0],
                   transfer_struct_c.ts_contribs_ss__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_intensity_db_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_intensity_db_,
      blitz::shape(transfer_struct_c.ts_intensity_db__f_shapes[0],
                   transfer_struct_c.ts_intensity_db__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_intensity_ss_.reference(blitz::Array<double, 3>(transfer_struct_c.ts_intensity_ss_,
      blitz::shape(transfer_struct_c.ts_intensity_ss__f_shapes[0],
                   transfer_struct_c.ts_intensity_ss__f_shapes[1],
                   transfer_struct_c.ts_intensity_ss__f_shapes[2]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<3>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_sup_ss transfer_struct_c;

  blitz::Array<double, 2> ts_contribs_ss_;
  blitz::Array<double, 2> ts_intensity_db_;
  blitz::Array<double, 3> ts_intensity_ss_;
  
};

// Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def_m" in file: "lidort_sup_def.f90"
extern "C" {
  void lidort_sup_inout_c_alloc_init(struct lidort_sup_inout *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_inout_c_init_only(struct lidort_sup_inout *transfer_struct_c, void **fortran_type_c);
  void lidort_sup_inout_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_sup_inout_c_destroy(void **fortran_type_c);
  
}

struct lidort_sup_inout {
  void* brdf_;
  int brdf__f_byte_size;

  void* sleave_;
  int sleave__f_byte_size;

  void* ss_;
  int ss__f_byte_size;

  
};

// Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def_m" in file: "lidort_sup_def.f90"
class Lidort_Sup_Inout : public Spurr_Sup_Inout_Base {
public:
  // Allocating constructor
  Lidort_Sup_Inout() : Spurr_Sup_Inout_Base() {
    lidort_sup_inout_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Sup_Inout(void* allocated_f_type_c) : Spurr_Sup_Inout_Base(allocated_f_type_c) {
    lidort_sup_inout_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Sup_Inout() {
    if (owns_pointer)
      lidort_sup_inout_c_destroy(&fortran_type_c);
  }

  Lidort_Sup_Brdf& brdf() {
    return *brdf_;
  }

  const Lidort_Sup_Brdf& brdf() const {
    return *brdf_;
  }

  void brdf(Lidort_Sup_Brdf& brdf_in) {
    void* src_ptr = brdf_in.fortran_type_ptr();
    void* dst_ptr = brdf_->fortran_type_ptr();
    lidort_sup_brdf_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Sup_Sleave& sleave() {
    return *sleave_;
  }

  const Lidort_Sup_Sleave& sleave() const {
    return *sleave_;
  }

  void sleave(Lidort_Sup_Sleave& sleave_in) {
    void* src_ptr = sleave_in.fortran_type_ptr();
    void* dst_ptr = sleave_->fortran_type_ptr();
    lidort_sup_sleave_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Sup_Ss& ss() {
    return *ss_;
  }

  const Lidort_Sup_Ss& ss() const {
    return *ss_;
  }

  void ss(Lidort_Sup_Ss& ss_in) {
    void* src_ptr = ss_in.fortran_type_ptr();
    void* dst_ptr = ss_->fortran_type_ptr();
    lidort_sup_ss_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Sup_Inout:" << std::endl
      << "  brdf: " << brdf()  << std::endl
      << "sleave: " << sleave()  << std::endl
      << "    ss: " << ss()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    brdf_.reset( new Lidort_Sup_Brdf(transfer_struct_c.brdf_) );
    sleave_.reset( new Lidort_Sup_Sleave(transfer_struct_c.sleave_) );
    ss_.reset( new Lidort_Sup_Ss(transfer_struct_c.ss_) );
    
  }

  struct lidort_sup_inout transfer_struct_c;

  boost::shared_ptr<Lidort_Sup_Brdf> brdf_;
  boost::shared_ptr<Lidort_Sup_Sleave> sleave_;
  boost::shared_ptr<Lidort_Sup_Ss> ss_;
  
};

// Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_fixed_boolean_c_alloc_init(struct lidort_fixed_boolean *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_boolean_c_init_only(struct lidort_fixed_boolean *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_boolean_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_boolean_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_boolean {
  int* ts_do_albtrn_media_;
  int ts_do_albtrn_media__f_shapes[1];
  int ts_do_albtrn_media__f_byte_size;

  int* ts_do_boa_illumination_;
  int ts_do_boa_illumination__f_byte_size;

  int* ts_do_brdf_surface_;
  int ts_do_brdf_surface__f_byte_size;

  int* ts_do_dnwelling_;
  int ts_do_dnwelling__f_byte_size;

  int* ts_do_fluorescence_;
  int ts_do_fluorescence__f_byte_size;

  int* ts_do_fullrad_mode_;
  int ts_do_fullrad_mode__f_byte_size;

  int* ts_do_mssts_;
  int ts_do_mssts__f_byte_size;

  int* ts_do_plane_parallel_;
  int ts_do_plane_parallel__f_byte_size;

  int* ts_do_planetary_problem_;
  int ts_do_planetary_problem__f_byte_size;

  int* ts_do_sl_isotropic_;
  int ts_do_sl_isotropic__f_byte_size;

  int* ts_do_surface_emission_;
  int ts_do_surface_emission__f_byte_size;

  int* ts_do_surface_leaving_;
  int ts_do_surface_leaving__f_byte_size;

  int* ts_do_tf_iteration_;
  int ts_do_tf_iteration__f_byte_size;

  int* ts_do_thermal_emission_;
  int ts_do_thermal_emission__f_byte_size;

  int* ts_do_toa_contribs_;
  int ts_do_toa_contribs__f_byte_size;

  int* ts_do_toa_illumination_;
  int ts_do_toa_illumination__f_byte_size;

  int* ts_do_upwelling_;
  int ts_do_upwelling__f_byte_size;

  int* ts_do_water_leaving_;
  int ts_do_water_leaving__f_byte_size;

  int* ts_do_wladjusted_output_;
  int ts_do_wladjusted_output__f_byte_size;

  
};

// Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Fixed_Boolean : public Spurr_Fixed_Boolean_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Boolean() : Spurr_Fixed_Boolean_Base() {
    lidort_fixed_boolean_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Boolean(void* allocated_f_type_c) : Spurr_Fixed_Boolean_Base(allocated_f_type_c) {
    lidort_fixed_boolean_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Boolean() {
    if (owns_pointer)
      lidort_fixed_boolean_c_destroy(&fortran_type_c);
  }

  const blitz::Array<bool, 1> ts_do_albtrn_media() const {
    blitz::Array<bool,1> as_bool(ts_do_albtrn_media_.shape());
    as_bool = blitz::where(ts_do_albtrn_media_ != 0, true, false);
    return as_bool;
  }

  void ts_do_albtrn_media(const blitz::Array<bool, 1>& ts_do_albtrn_media_in) {
    blitz::Array<int,1> as_int(ts_do_albtrn_media_.shape());
    as_int = blitz::where(ts_do_albtrn_media_in == true, FORTRAN_TRUE_INT, 0);
    ts_do_albtrn_media_ = as_int;
  }

  
  const bool ts_do_boa_illumination() const {
    return *transfer_struct_c.ts_do_boa_illumination_ != 0;
  }

  void ts_do_boa_illumination(const bool& ts_do_boa_illumination_in) {
    *transfer_struct_c.ts_do_boa_illumination_ = ts_do_boa_illumination_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_brdf_surface() const {
    return *transfer_struct_c.ts_do_brdf_surface_ != 0;
  }

  void ts_do_brdf_surface(const bool& ts_do_brdf_surface_in) {
    *transfer_struct_c.ts_do_brdf_surface_ = ts_do_brdf_surface_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_dnwelling() const {
    return *transfer_struct_c.ts_do_dnwelling_ != 0;
  }

  void ts_do_dnwelling(const bool& ts_do_dnwelling_in) {
    *transfer_struct_c.ts_do_dnwelling_ = ts_do_dnwelling_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_fluorescence() const {
    return *transfer_struct_c.ts_do_fluorescence_ != 0;
  }

  void ts_do_fluorescence(const bool& ts_do_fluorescence_in) {
    *transfer_struct_c.ts_do_fluorescence_ = ts_do_fluorescence_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_fullrad_mode() const {
    return *transfer_struct_c.ts_do_fullrad_mode_ != 0;
  }

  void ts_do_fullrad_mode(const bool& ts_do_fullrad_mode_in) {
    *transfer_struct_c.ts_do_fullrad_mode_ = ts_do_fullrad_mode_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_mssts() const {
    return *transfer_struct_c.ts_do_mssts_ != 0;
  }

  void ts_do_mssts(const bool& ts_do_mssts_in) {
    *transfer_struct_c.ts_do_mssts_ = ts_do_mssts_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_plane_parallel() const {
    return *transfer_struct_c.ts_do_plane_parallel_ != 0;
  }

  void ts_do_plane_parallel(const bool& ts_do_plane_parallel_in) {
    *transfer_struct_c.ts_do_plane_parallel_ = ts_do_plane_parallel_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_planetary_problem() const {
    return *transfer_struct_c.ts_do_planetary_problem_ != 0;
  }

  void ts_do_planetary_problem(const bool& ts_do_planetary_problem_in) {
    *transfer_struct_c.ts_do_planetary_problem_ = ts_do_planetary_problem_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_sl_isotropic() const {
    return *transfer_struct_c.ts_do_sl_isotropic_ != 0;
  }

  void ts_do_sl_isotropic(const bool& ts_do_sl_isotropic_in) {
    *transfer_struct_c.ts_do_sl_isotropic_ = ts_do_sl_isotropic_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_surface_emission() const {
    return *transfer_struct_c.ts_do_surface_emission_ != 0;
  }

  void ts_do_surface_emission(const bool& ts_do_surface_emission_in) {
    *transfer_struct_c.ts_do_surface_emission_ = ts_do_surface_emission_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_surface_leaving() const {
    return *transfer_struct_c.ts_do_surface_leaving_ != 0;
  }

  void ts_do_surface_leaving(const bool& ts_do_surface_leaving_in) {
    *transfer_struct_c.ts_do_surface_leaving_ = ts_do_surface_leaving_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_tf_iteration() const {
    return *transfer_struct_c.ts_do_tf_iteration_ != 0;
  }

  void ts_do_tf_iteration(const bool& ts_do_tf_iteration_in) {
    *transfer_struct_c.ts_do_tf_iteration_ = ts_do_tf_iteration_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_thermal_emission() const {
    return *transfer_struct_c.ts_do_thermal_emission_ != 0;
  }

  void ts_do_thermal_emission(const bool& ts_do_thermal_emission_in) {
    *transfer_struct_c.ts_do_thermal_emission_ = ts_do_thermal_emission_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_toa_contribs() const {
    return *transfer_struct_c.ts_do_toa_contribs_ != 0;
  }

  void ts_do_toa_contribs(const bool& ts_do_toa_contribs_in) {
    *transfer_struct_c.ts_do_toa_contribs_ = ts_do_toa_contribs_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_toa_illumination() const {
    return *transfer_struct_c.ts_do_toa_illumination_ != 0;
  }

  void ts_do_toa_illumination(const bool& ts_do_toa_illumination_in) {
    *transfer_struct_c.ts_do_toa_illumination_ = ts_do_toa_illumination_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_upwelling() const {
    return *transfer_struct_c.ts_do_upwelling_ != 0;
  }

  void ts_do_upwelling(const bool& ts_do_upwelling_in) {
    *transfer_struct_c.ts_do_upwelling_ = ts_do_upwelling_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_water_leaving() const {
    return *transfer_struct_c.ts_do_water_leaving_ != 0;
  }

  void ts_do_water_leaving(const bool& ts_do_water_leaving_in) {
    *transfer_struct_c.ts_do_water_leaving_ = ts_do_water_leaving_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_wladjusted_output() const {
    return *transfer_struct_c.ts_do_wladjusted_output_ != 0;
  }

  void ts_do_wladjusted_output(const bool& ts_do_wladjusted_output_in) {
    *transfer_struct_c.ts_do_wladjusted_output_ = ts_do_wladjusted_output_in ? FORTRAN_TRUE_INT : 0;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Boolean:" << std::endl
      << "     ts_do_albtrn_media: " << std::endl << ts_do_albtrn_media()  << std::endl
      << " ts_do_boa_illumination: " << ts_do_boa_illumination()  << std::endl
      << "     ts_do_brdf_surface: " << ts_do_brdf_surface()  << std::endl
      << "        ts_do_dnwelling: " << ts_do_dnwelling()  << std::endl
      << "     ts_do_fluorescence: " << ts_do_fluorescence()  << std::endl
      << "     ts_do_fullrad_mode: " << ts_do_fullrad_mode()  << std::endl
      << "            ts_do_mssts: " << ts_do_mssts()  << std::endl
      << "   ts_do_plane_parallel: " << ts_do_plane_parallel()  << std::endl
      << "ts_do_planetary_problem: " << ts_do_planetary_problem()  << std::endl
      << "     ts_do_sl_isotropic: " << ts_do_sl_isotropic()  << std::endl
      << " ts_do_surface_emission: " << ts_do_surface_emission()  << std::endl
      << "  ts_do_surface_leaving: " << ts_do_surface_leaving()  << std::endl
      << "     ts_do_tf_iteration: " << ts_do_tf_iteration()  << std::endl
      << " ts_do_thermal_emission: " << ts_do_thermal_emission()  << std::endl
      << "     ts_do_toa_contribs: " << ts_do_toa_contribs()  << std::endl
      << " ts_do_toa_illumination: " << ts_do_toa_illumination()  << std::endl
      << "        ts_do_upwelling: " << ts_do_upwelling()  << std::endl
      << "    ts_do_water_leaving: " << ts_do_water_leaving()  << std::endl
      << "ts_do_wladjusted_output: " << ts_do_wladjusted_output()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_do_albtrn_media_",sizeof(*transfer_struct_c.ts_do_albtrn_media_),transfer_struct_c.ts_do_albtrn_media__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_boa_illumination_",sizeof(*transfer_struct_c.ts_do_boa_illumination_),transfer_struct_c.ts_do_boa_illumination__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_brdf_surface_",sizeof(*transfer_struct_c.ts_do_brdf_surface_),transfer_struct_c.ts_do_brdf_surface__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_dnwelling_",sizeof(*transfer_struct_c.ts_do_dnwelling_),transfer_struct_c.ts_do_dnwelling__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_fluorescence_",sizeof(*transfer_struct_c.ts_do_fluorescence_),transfer_struct_c.ts_do_fluorescence__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_fullrad_mode_",sizeof(*transfer_struct_c.ts_do_fullrad_mode_),transfer_struct_c.ts_do_fullrad_mode__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_mssts_",sizeof(*transfer_struct_c.ts_do_mssts_),transfer_struct_c.ts_do_mssts__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_plane_parallel_",sizeof(*transfer_struct_c.ts_do_plane_parallel_),transfer_struct_c.ts_do_plane_parallel__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_planetary_problem_",sizeof(*transfer_struct_c.ts_do_planetary_problem_),transfer_struct_c.ts_do_planetary_problem__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_sl_isotropic_",sizeof(*transfer_struct_c.ts_do_sl_isotropic_),transfer_struct_c.ts_do_sl_isotropic__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_surface_emission_",sizeof(*transfer_struct_c.ts_do_surface_emission_),transfer_struct_c.ts_do_surface_emission__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_surface_leaving_",sizeof(*transfer_struct_c.ts_do_surface_leaving_),transfer_struct_c.ts_do_surface_leaving__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_tf_iteration_",sizeof(*transfer_struct_c.ts_do_tf_iteration_),transfer_struct_c.ts_do_tf_iteration__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_thermal_emission_",sizeof(*transfer_struct_c.ts_do_thermal_emission_),transfer_struct_c.ts_do_thermal_emission__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_toa_contribs_",sizeof(*transfer_struct_c.ts_do_toa_contribs_),transfer_struct_c.ts_do_toa_contribs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_toa_illumination_",sizeof(*transfer_struct_c.ts_do_toa_illumination_),transfer_struct_c.ts_do_toa_illumination__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_upwelling_",sizeof(*transfer_struct_c.ts_do_upwelling_),transfer_struct_c.ts_do_upwelling__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_water_leaving_",sizeof(*transfer_struct_c.ts_do_water_leaving_),transfer_struct_c.ts_do_water_leaving__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_wladjusted_output_",sizeof(*transfer_struct_c.ts_do_wladjusted_output_),transfer_struct_c.ts_do_wladjusted_output__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_do_albtrn_media_.reference(blitz::Array<int, 1>(transfer_struct_c.ts_do_albtrn_media_,
      blitz::shape(transfer_struct_c.ts_do_albtrn_media__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_boolean transfer_struct_c;

  blitz::Array<int, 1> ts_do_albtrn_media_;
  
};

// Links to type: "lidort_fixed_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_fixed_control_c_alloc_init(struct lidort_fixed_control *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_control_c_init_only(struct lidort_fixed_control *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_control_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_control_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_control {
  double* ts_asymtx_tolerance_;
  int ts_asymtx_tolerance__f_byte_size;

  double* ts_boa_illumination_;
  int ts_boa_illumination__f_byte_size;

  double* ts_lidort_accuracy_;
  int ts_lidort_accuracy__f_byte_size;

  int* ts_n_thermal_coeffs_;
  int ts_n_thermal_coeffs__f_byte_size;

  int* ts_nfinelayers_;
  int ts_nfinelayers__f_byte_size;

  int* ts_nlayers_;
  int ts_nlayers__f_byte_size;

  int* ts_nstreams_;
  int ts_nstreams__f_byte_size;

  int* ts_taylor_order_;
  int ts_taylor_order__f_byte_size;

  double* ts_tf_criterion_;
  int ts_tf_criterion__f_byte_size;

  int* ts_tf_maxiter_;
  int ts_tf_maxiter__f_byte_size;

  double* ts_toa_illumination_;
  int ts_toa_illumination__f_byte_size;

  
};

// Links to type: "lidort_fixed_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Fixed_Control : public Spurr_Fixed_Control_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Control() : Spurr_Fixed_Control_Base() {
    lidort_fixed_control_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Control(void* allocated_f_type_c) : Spurr_Fixed_Control_Base(allocated_f_type_c) {
    lidort_fixed_control_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Control() {
    if (owns_pointer)
      lidort_fixed_control_c_destroy(&fortran_type_c);
  }

  const double& ts_asymtx_tolerance() const {
    return *transfer_struct_c.ts_asymtx_tolerance_;
  }

  void ts_asymtx_tolerance(const double& ts_asymtx_tolerance_in) {
    *transfer_struct_c.ts_asymtx_tolerance_ = ts_asymtx_tolerance_in;
  }

  
  const double& ts_boa_illumination() const {
    return *transfer_struct_c.ts_boa_illumination_;
  }

  void ts_boa_illumination(const double& ts_boa_illumination_in) {
    *transfer_struct_c.ts_boa_illumination_ = ts_boa_illumination_in;
  }

  
  const double& ts_lidort_accuracy() const {
    return *transfer_struct_c.ts_lidort_accuracy_;
  }

  void ts_lidort_accuracy(const double& ts_lidort_accuracy_in) {
    *transfer_struct_c.ts_lidort_accuracy_ = ts_lidort_accuracy_in;
  }

  
  const int& ts_n_thermal_coeffs() const {
    return *transfer_struct_c.ts_n_thermal_coeffs_;
  }

  void ts_n_thermal_coeffs(const int& ts_n_thermal_coeffs_in) {
    *transfer_struct_c.ts_n_thermal_coeffs_ = ts_n_thermal_coeffs_in;
  }

  
  const int& ts_nfinelayers() const {
    return *transfer_struct_c.ts_nfinelayers_;
  }

  void ts_nfinelayers(const int& ts_nfinelayers_in) {
    *transfer_struct_c.ts_nfinelayers_ = ts_nfinelayers_in;
  }

  
  const int& ts_nlayers() const {
    return *transfer_struct_c.ts_nlayers_;
  }

  void ts_nlayers(const int& ts_nlayers_in) {
    *transfer_struct_c.ts_nlayers_ = ts_nlayers_in;
  }

  
  const int& ts_nstreams() const {
    return *transfer_struct_c.ts_nstreams_;
  }

  void ts_nstreams(const int& ts_nstreams_in) {
    *transfer_struct_c.ts_nstreams_ = ts_nstreams_in;
  }

  
  const int& ts_taylor_order() const {
    return *transfer_struct_c.ts_taylor_order_;
  }

  void ts_taylor_order(const int& ts_taylor_order_in) {
    *transfer_struct_c.ts_taylor_order_ = ts_taylor_order_in;
  }

  
  const double& ts_tf_criterion() const {
    return *transfer_struct_c.ts_tf_criterion_;
  }

  void ts_tf_criterion(const double& ts_tf_criterion_in) {
    *transfer_struct_c.ts_tf_criterion_ = ts_tf_criterion_in;
  }

  
  const int& ts_tf_maxiter() const {
    return *transfer_struct_c.ts_tf_maxiter_;
  }

  void ts_tf_maxiter(const int& ts_tf_maxiter_in) {
    *transfer_struct_c.ts_tf_maxiter_ = ts_tf_maxiter_in;
  }

  
  const double& ts_toa_illumination() const {
    return *transfer_struct_c.ts_toa_illumination_;
  }

  void ts_toa_illumination(const double& ts_toa_illumination_in) {
    *transfer_struct_c.ts_toa_illumination_ = ts_toa_illumination_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Control:" << std::endl
      << "ts_asymtx_tolerance: " << ts_asymtx_tolerance()  << std::endl
      << "ts_boa_illumination: " << ts_boa_illumination()  << std::endl
      << " ts_lidort_accuracy: " << ts_lidort_accuracy()  << std::endl
      << "ts_n_thermal_coeffs: " << ts_n_thermal_coeffs()  << std::endl
      << "     ts_nfinelayers: " << ts_nfinelayers()  << std::endl
      << "         ts_nlayers: " << ts_nlayers()  << std::endl
      << "        ts_nstreams: " << ts_nstreams()  << std::endl
      << "    ts_taylor_order: " << ts_taylor_order()  << std::endl
      << "    ts_tf_criterion: " << ts_tf_criterion()  << std::endl
      << "      ts_tf_maxiter: " << ts_tf_maxiter()  << std::endl
      << "ts_toa_illumination: " << ts_toa_illumination()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_asymtx_tolerance_",sizeof(*transfer_struct_c.ts_asymtx_tolerance_),transfer_struct_c.ts_asymtx_tolerance__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_boa_illumination_",sizeof(*transfer_struct_c.ts_boa_illumination_),transfer_struct_c.ts_boa_illumination__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lidort_accuracy_",sizeof(*transfer_struct_c.ts_lidort_accuracy_),transfer_struct_c.ts_lidort_accuracy__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_thermal_coeffs_",sizeof(*transfer_struct_c.ts_n_thermal_coeffs_),transfer_struct_c.ts_n_thermal_coeffs__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_nfinelayers_",sizeof(*transfer_struct_c.ts_nfinelayers_),transfer_struct_c.ts_nfinelayers__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_nlayers_",sizeof(*transfer_struct_c.ts_nlayers_),transfer_struct_c.ts_nlayers__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_nstreams_",sizeof(*transfer_struct_c.ts_nstreams_),transfer_struct_c.ts_nstreams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_taylor_order_",sizeof(*transfer_struct_c.ts_taylor_order_),transfer_struct_c.ts_taylor_order__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_tf_criterion_",sizeof(*transfer_struct_c.ts_tf_criterion_),transfer_struct_c.ts_tf_criterion__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_tf_maxiter_",sizeof(*transfer_struct_c.ts_tf_maxiter_),transfer_struct_c.ts_tf_maxiter__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_toa_illumination_",sizeof(*transfer_struct_c.ts_toa_illumination_),transfer_struct_c.ts_toa_illumination__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_control transfer_struct_c;

  
};

// Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_fixed_sunrays_c_alloc_init(struct lidort_fixed_sunrays *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_sunrays_c_init_only(struct lidort_fixed_sunrays *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_sunrays_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_sunrays_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_sunrays {
  double* ts_flux_factor_;
  int ts_flux_factor__f_byte_size;

  
};

// Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Fixed_Sunrays : public Spurr_Fixed_Sunrays_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Sunrays() : Spurr_Fixed_Sunrays_Base() {
    lidort_fixed_sunrays_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Sunrays(void* allocated_f_type_c) : Spurr_Fixed_Sunrays_Base(allocated_f_type_c) {
    lidort_fixed_sunrays_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Sunrays() {
    if (owns_pointer)
      lidort_fixed_sunrays_c_destroy(&fortran_type_c);
  }

  const double& ts_flux_factor() const {
    return *transfer_struct_c.ts_flux_factor_;
  }

  void ts_flux_factor(const double& ts_flux_factor_in) {
    *transfer_struct_c.ts_flux_factor_ = ts_flux_factor_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Sunrays:" << std::endl
      << "ts_flux_factor: " << ts_flux_factor()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_flux_factor_",sizeof(*transfer_struct_c.ts_flux_factor_),transfer_struct_c.ts_flux_factor__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_sunrays transfer_struct_c;

  
};

// Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_fixed_uservalues_c_alloc_init(struct lidort_fixed_uservalues *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_uservalues_c_init_only(struct lidort_fixed_uservalues *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_uservalues_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_uservalues_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_uservalues {
  int* ts_n_user_levels_;
  int ts_n_user_levels__f_byte_size;

  
};

// Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Fixed_Uservalues : public Spurr_Fixed_Uservalues_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Uservalues() : Spurr_Fixed_Uservalues_Base() {
    lidort_fixed_uservalues_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Uservalues(void* allocated_f_type_c) : Spurr_Fixed_Uservalues_Base(allocated_f_type_c) {
    lidort_fixed_uservalues_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Uservalues() {
    if (owns_pointer)
      lidort_fixed_uservalues_c_destroy(&fortran_type_c);
  }

  const int& ts_n_user_levels() const {
    return *transfer_struct_c.ts_n_user_levels_;
  }

  void ts_n_user_levels(const int& ts_n_user_levels_in) {
    *transfer_struct_c.ts_n_user_levels_ = ts_n_user_levels_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Uservalues:" << std::endl
      << "ts_n_user_levels: " << ts_n_user_levels()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_n_user_levels_",sizeof(*transfer_struct_c.ts_n_user_levels_),transfer_struct_c.ts_n_user_levels__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_uservalues transfer_struct_c;

  
};

// Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_fixed_chapman_c_alloc_init(struct lidort_fixed_chapman *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_chapman_c_init_only(struct lidort_fixed_chapman *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_chapman_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_chapman_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_chapman {
  int* ts_finegrid_;
  int ts_finegrid__f_shapes[1];
  int ts_finegrid__f_byte_size;

  double* ts_height_grid_;
  int ts_height_grid__f_shapes[1];
  int ts_height_grid__f_byte_size;

  double* ts_pressure_grid_;
  int ts_pressure_grid__f_shapes[1];
  int ts_pressure_grid__f_byte_size;

  double* ts_rfindex_parameter_;
  int ts_rfindex_parameter__f_byte_size;

  double* ts_temperature_grid_;
  int ts_temperature_grid__f_shapes[1];
  int ts_temperature_grid__f_byte_size;

  
};

// Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Fixed_Chapman : public Spurr_Fixed_Chapman_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Chapman() : Spurr_Fixed_Chapman_Base() {
    lidort_fixed_chapman_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Chapman(void* allocated_f_type_c) : Spurr_Fixed_Chapman_Base(allocated_f_type_c) {
    lidort_fixed_chapman_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Chapman() {
    if (owns_pointer)
      lidort_fixed_chapman_c_destroy(&fortran_type_c);
  }

  const blitz::Array<int, 1>& ts_finegrid() const {
    return ts_finegrid_;
  }

  void ts_finegrid(const blitz::Array<int, 1>& ts_finegrid_in) {
    ts_finegrid_ = ts_finegrid_in;
  }

  
  const blitz::Array<double, 1>& ts_height_grid() const {
    return ts_height_grid_;
  }

  void ts_height_grid(const blitz::Array<double, 1>& ts_height_grid_in) {
    ts_height_grid_ = ts_height_grid_in;
  }

  
  const blitz::Array<double, 1>& ts_pressure_grid() const {
    return ts_pressure_grid_;
  }

  void ts_pressure_grid(const blitz::Array<double, 1>& ts_pressure_grid_in) {
    ts_pressure_grid_ = ts_pressure_grid_in;
  }

  
  const double& ts_rfindex_parameter() const {
    return *transfer_struct_c.ts_rfindex_parameter_;
  }

  void ts_rfindex_parameter(const double& ts_rfindex_parameter_in) {
    *transfer_struct_c.ts_rfindex_parameter_ = ts_rfindex_parameter_in;
  }

  
  const blitz::Array<double, 1>& ts_temperature_grid() const {
    return ts_temperature_grid_;
  }

  void ts_temperature_grid(const blitz::Array<double, 1>& ts_temperature_grid_in) {
    ts_temperature_grid_ = ts_temperature_grid_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Chapman:" << std::endl
      << "         ts_finegrid: " << std::endl << ts_finegrid()  << std::endl
      << "      ts_height_grid: " << std::endl << ts_height_grid()  << std::endl
      << "    ts_pressure_grid: " << std::endl << ts_pressure_grid()  << std::endl
      << "ts_rfindex_parameter: " << ts_rfindex_parameter()  << std::endl
      << " ts_temperature_grid: " << std::endl << ts_temperature_grid()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_finegrid_",sizeof(*transfer_struct_c.ts_finegrid_),transfer_struct_c.ts_finegrid__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_height_grid_",sizeof(*transfer_struct_c.ts_height_grid_),transfer_struct_c.ts_height_grid__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_pressure_grid_",sizeof(*transfer_struct_c.ts_pressure_grid_),transfer_struct_c.ts_pressure_grid__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_rfindex_parameter_",sizeof(*transfer_struct_c.ts_rfindex_parameter_),transfer_struct_c.ts_rfindex_parameter__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_temperature_grid_",sizeof(*transfer_struct_c.ts_temperature_grid_),transfer_struct_c.ts_temperature_grid__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_finegrid_.reference(blitz::Array<int, 1>(transfer_struct_c.ts_finegrid_,
      blitz::shape(transfer_struct_c.ts_finegrid__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_height_grid_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_height_grid_,
      blitz::shape(transfer_struct_c.ts_height_grid__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_pressure_grid_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_pressure_grid_,
      blitz::shape(transfer_struct_c.ts_pressure_grid__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_temperature_grid_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_temperature_grid_,
      blitz::shape(transfer_struct_c.ts_temperature_grid__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_chapman transfer_struct_c;

  blitz::Array<int, 1> ts_finegrid_;
  blitz::Array<double, 1> ts_height_grid_;
  blitz::Array<double, 1> ts_pressure_grid_;
  blitz::Array<double, 1> ts_temperature_grid_;
  
};

// Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_fixed_optical_c_alloc_init(struct lidort_fixed_optical *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_optical_c_init_only(struct lidort_fixed_optical *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_optical_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_optical_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_optical {
  double* ts_atmos_wavelength_;
  int ts_atmos_wavelength__f_byte_size;

  double* ts_deltau_vert_input_;
  int ts_deltau_vert_input__f_shapes[1];
  int ts_deltau_vert_input__f_byte_size;

  double* ts_lambertian_albedo_;
  int ts_lambertian_albedo__f_byte_size;

  double* ts_phasfunc_input_dn_;
  int ts_phasfunc_input_dn__f_shapes[2];
  int ts_phasfunc_input_dn__f_byte_size;

  double* ts_phasfunc_input_up_;
  int ts_phasfunc_input_up__f_shapes[2];
  int ts_phasfunc_input_up__f_byte_size;

  double* ts_phasmoms_total_input_;
  int ts_phasmoms_total_input__f_shapes[2];
  int ts_phasmoms_total_input__f_byte_size;

  double* ts_surface_bb_input_;
  int ts_surface_bb_input__f_byte_size;

  double* ts_thermal_bb_input_;
  int ts_thermal_bb_input__f_shapes[1];
  int ts_thermal_bb_input__f_byte_size;

  
};

// Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Fixed_Optical : public Spurr_Fixed_Optical_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Optical() : Spurr_Fixed_Optical_Base() {
    lidort_fixed_optical_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Optical(void* allocated_f_type_c) : Spurr_Fixed_Optical_Base(allocated_f_type_c) {
    lidort_fixed_optical_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Optical() {
    if (owns_pointer)
      lidort_fixed_optical_c_destroy(&fortran_type_c);
  }

  const double& ts_atmos_wavelength() const {
    return *transfer_struct_c.ts_atmos_wavelength_;
  }

  void ts_atmos_wavelength(const double& ts_atmos_wavelength_in) {
    *transfer_struct_c.ts_atmos_wavelength_ = ts_atmos_wavelength_in;
  }

  
  const blitz::Array<double, 1>& ts_deltau_vert_input() const {
    return ts_deltau_vert_input_;
  }

  void ts_deltau_vert_input(const blitz::Array<double, 1>& ts_deltau_vert_input_in) {
    ts_deltau_vert_input_ = ts_deltau_vert_input_in;
  }

  
  const double& ts_lambertian_albedo() const {
    return *transfer_struct_c.ts_lambertian_albedo_;
  }

  void ts_lambertian_albedo(const double& ts_lambertian_albedo_in) {
    *transfer_struct_c.ts_lambertian_albedo_ = ts_lambertian_albedo_in;
  }

  
  const blitz::Array<double, 2>& ts_phasfunc_input_dn() const {
    return ts_phasfunc_input_dn_;
  }

  void ts_phasfunc_input_dn(const blitz::Array<double, 2>& ts_phasfunc_input_dn_in) {
    ts_phasfunc_input_dn_ = ts_phasfunc_input_dn_in;
  }

  
  const blitz::Array<double, 2>& ts_phasfunc_input_up() const {
    return ts_phasfunc_input_up_;
  }

  void ts_phasfunc_input_up(const blitz::Array<double, 2>& ts_phasfunc_input_up_in) {
    ts_phasfunc_input_up_ = ts_phasfunc_input_up_in;
  }

  
  const blitz::Array<double, 2>& ts_phasmoms_total_input() const {
    return ts_phasmoms_total_input_;
  }

  void ts_phasmoms_total_input(const blitz::Array<double, 2>& ts_phasmoms_total_input_in) {
    ts_phasmoms_total_input_ = ts_phasmoms_total_input_in;
  }

  
  const double& ts_surface_bb_input() const {
    return *transfer_struct_c.ts_surface_bb_input_;
  }

  void ts_surface_bb_input(const double& ts_surface_bb_input_in) {
    *transfer_struct_c.ts_surface_bb_input_ = ts_surface_bb_input_in;
  }

  
  const blitz::Array<double, 1>& ts_thermal_bb_input() const {
    return ts_thermal_bb_input_;
  }

  void ts_thermal_bb_input(const blitz::Array<double, 1>& ts_thermal_bb_input_in) {
    ts_thermal_bb_input_ = ts_thermal_bb_input_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Optical:" << std::endl
      << "    ts_atmos_wavelength: " << ts_atmos_wavelength()  << std::endl
      << "   ts_deltau_vert_input: " << std::endl << ts_deltau_vert_input()  << std::endl
      << "   ts_lambertian_albedo: " << ts_lambertian_albedo()  << std::endl
      << "   ts_phasfunc_input_dn: " << std::endl << ts_phasfunc_input_dn()  << std::endl
      << "   ts_phasfunc_input_up: " << std::endl << ts_phasfunc_input_up()  << std::endl
      << "ts_phasmoms_total_input: " << std::endl << ts_phasmoms_total_input()  << std::endl
      << "    ts_surface_bb_input: " << ts_surface_bb_input()  << std::endl
      << "    ts_thermal_bb_input: " << std::endl << ts_thermal_bb_input()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_atmos_wavelength_",sizeof(*transfer_struct_c.ts_atmos_wavelength_),transfer_struct_c.ts_atmos_wavelength__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_deltau_vert_input_",sizeof(*transfer_struct_c.ts_deltau_vert_input_),transfer_struct_c.ts_deltau_vert_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_lambertian_albedo_",sizeof(*transfer_struct_c.ts_lambertian_albedo_),transfer_struct_c.ts_lambertian_albedo__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_phasfunc_input_dn_",sizeof(*transfer_struct_c.ts_phasfunc_input_dn_),transfer_struct_c.ts_phasfunc_input_dn__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_phasfunc_input_up_",sizeof(*transfer_struct_c.ts_phasfunc_input_up_),transfer_struct_c.ts_phasfunc_input_up__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_phasmoms_total_input_",sizeof(*transfer_struct_c.ts_phasmoms_total_input_),transfer_struct_c.ts_phasmoms_total_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_surface_bb_input_",sizeof(*transfer_struct_c.ts_surface_bb_input_),transfer_struct_c.ts_surface_bb_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_thermal_bb_input_",sizeof(*transfer_struct_c.ts_thermal_bb_input_),transfer_struct_c.ts_thermal_bb_input__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_deltau_vert_input_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_deltau_vert_input_,
      blitz::shape(transfer_struct_c.ts_deltau_vert_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_phasfunc_input_dn_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_phasfunc_input_dn_,
      blitz::shape(transfer_struct_c.ts_phasfunc_input_dn__f_shapes[0],
                   transfer_struct_c.ts_phasfunc_input_dn__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_phasfunc_input_up_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_phasfunc_input_up_,
      blitz::shape(transfer_struct_c.ts_phasfunc_input_up__f_shapes[0],
                   transfer_struct_c.ts_phasfunc_input_up__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_phasmoms_total_input_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_phasmoms_total_input_,
      blitz::shape(transfer_struct_c.ts_phasmoms_total_input__f_shapes[0],
                   transfer_struct_c.ts_phasmoms_total_input__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_thermal_bb_input_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_thermal_bb_input_,
      blitz::shape(transfer_struct_c.ts_thermal_bb_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_optical transfer_struct_c;

  blitz::Array<double, 1> ts_deltau_vert_input_;
  blitz::Array<double, 2> ts_phasfunc_input_dn_;
  blitz::Array<double, 2> ts_phasfunc_input_up_;
  blitz::Array<double, 2> ts_phasmoms_total_input_;
  blitz::Array<double, 1> ts_thermal_bb_input_;
  
};

// Links to type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_fixed_write_c_alloc_init(struct lidort_fixed_write *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_write_c_init_only(struct lidort_fixed_write *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_write_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_write_c_destroy(void **fortran_type_c);
  void fixed_write_ts_fourier_write_filename_get(void **fortran_type_c, const int* ts_fourier_write_filename_in_len, const char* ts_fourier_write_filename_in);
  void fixed_write_ts_input_write_filename_get(void **fortran_type_c, const int* ts_input_write_filename_in_len, const char* ts_input_write_filename_in);
  void fixed_write_ts_results_write_filename_get(void **fortran_type_c, const int* ts_results_write_filename_in_len, const char* ts_results_write_filename_in);
  void fixed_write_ts_scenario_write_filename_get(void **fortran_type_c, const int* ts_scenario_write_filename_in_len, const char* ts_scenario_write_filename_in);
  
}

struct lidort_fixed_write {
  int* ts_do_debug_write_;
  int ts_do_debug_write__f_byte_size;

  int* ts_do_write_fourier_;
  int ts_do_write_fourier__f_byte_size;

  int* ts_do_write_input_;
  int ts_do_write_input__f_byte_size;

  int* ts_do_write_results_;
  int ts_do_write_results__f_byte_size;

  int* ts_do_write_scenario_;
  int ts_do_write_scenario__f_byte_size;

  
  int ts_fourier_write_filename__f_len;

  
  int ts_input_write_filename__f_len;

  
  int ts_results_write_filename__f_len;

  
  int ts_scenario_write_filename__f_len;

  
};

// Links to type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Fixed_Write : public Spurr_Fixed_Write_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Write() : Spurr_Fixed_Write_Base() {
    lidort_fixed_write_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Write(void* allocated_f_type_c) : Spurr_Fixed_Write_Base(allocated_f_type_c) {
    lidort_fixed_write_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Write() {
    if (owns_pointer)
      lidort_fixed_write_c_destroy(&fortran_type_c);
  }

  const bool ts_do_debug_write() const {
    return *transfer_struct_c.ts_do_debug_write_ != 0;
  }

  void ts_do_debug_write(const bool& ts_do_debug_write_in) {
    *transfer_struct_c.ts_do_debug_write_ = ts_do_debug_write_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_write_fourier() const {
    return *transfer_struct_c.ts_do_write_fourier_ != 0;
  }

  void ts_do_write_fourier(const bool& ts_do_write_fourier_in) {
    *transfer_struct_c.ts_do_write_fourier_ = ts_do_write_fourier_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_write_input() const {
    return *transfer_struct_c.ts_do_write_input_ != 0;
  }

  void ts_do_write_input(const bool& ts_do_write_input_in) {
    *transfer_struct_c.ts_do_write_input_ = ts_do_write_input_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_write_results() const {
    return *transfer_struct_c.ts_do_write_results_ != 0;
  }

  void ts_do_write_results(const bool& ts_do_write_results_in) {
    *transfer_struct_c.ts_do_write_results_ = ts_do_write_results_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_write_scenario() const {
    return *transfer_struct_c.ts_do_write_scenario_ != 0;
  }

  void ts_do_write_scenario(const bool& ts_do_write_scenario_in) {
    *transfer_struct_c.ts_do_write_scenario_ = ts_do_write_scenario_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const std::string ts_fourier_write_filename() const {
    std::string ts_fourier_write_filename_ret;
    blitz::Array<char, 1> ts_fourier_write_filename_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_fourier_write_filename__f_len+1, blitz::ColumnMajorArray<1>());
    fixed_write_ts_fourier_write_filename_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_fourier_write_filename__f_len, ts_fourier_write_filename_lcl.dataFirst());
    ts_fourier_write_filename_ret = ( std::string(std::string(ts_fourier_write_filename_lcl(blitz::Range::all()).begin(), ts_fourier_write_filename_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_fourier_write_filename_ret;
  }

  
  const std::string ts_input_write_filename() const {
    std::string ts_input_write_filename_ret;
    blitz::Array<char, 1> ts_input_write_filename_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_input_write_filename__f_len+1, blitz::ColumnMajorArray<1>());
    fixed_write_ts_input_write_filename_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_input_write_filename__f_len, ts_input_write_filename_lcl.dataFirst());
    ts_input_write_filename_ret = ( std::string(std::string(ts_input_write_filename_lcl(blitz::Range::all()).begin(), ts_input_write_filename_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_input_write_filename_ret;
  }

  
  const std::string ts_results_write_filename() const {
    std::string ts_results_write_filename_ret;
    blitz::Array<char, 1> ts_results_write_filename_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_results_write_filename__f_len+1, blitz::ColumnMajorArray<1>());
    fixed_write_ts_results_write_filename_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_results_write_filename__f_len, ts_results_write_filename_lcl.dataFirst());
    ts_results_write_filename_ret = ( std::string(std::string(ts_results_write_filename_lcl(blitz::Range::all()).begin(), ts_results_write_filename_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_results_write_filename_ret;
  }

  
  const std::string ts_scenario_write_filename() const {
    std::string ts_scenario_write_filename_ret;
    blitz::Array<char, 1> ts_scenario_write_filename_lcl = blitz::Array<char, 1>(transfer_struct_c.ts_scenario_write_filename__f_len+1, blitz::ColumnMajorArray<1>());
    fixed_write_ts_scenario_write_filename_get(const_cast<void**>(&fortran_type_c), &transfer_struct_c.ts_scenario_write_filename__f_len, ts_scenario_write_filename_lcl.dataFirst());
    ts_scenario_write_filename_ret = ( std::string(std::string(ts_scenario_write_filename_lcl(blitz::Range::all()).begin(), ts_scenario_write_filename_lcl(blitz::Range::all()).end()).c_str()) );
    return ts_scenario_write_filename_ret;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Write:" << std::endl
      << "         ts_do_debug_write: " << ts_do_debug_write()  << std::endl
      << "       ts_do_write_fourier: " << ts_do_write_fourier()  << std::endl
      << "         ts_do_write_input: " << ts_do_write_input()  << std::endl
      << "       ts_do_write_results: " << ts_do_write_results()  << std::endl
      << "      ts_do_write_scenario: " << ts_do_write_scenario()  << std::endl
      << " ts_fourier_write_filename: " << "\"" << ts_fourier_write_filename() << "\"" << std::endl
      << "   ts_input_write_filename: " << "\"" << ts_input_write_filename() << "\"" << std::endl
      << " ts_results_write_filename: " << "\"" << ts_results_write_filename() << "\"" << std::endl
      << "ts_scenario_write_filename: " << "\"" << ts_scenario_write_filename() << "\"" << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_do_debug_write_",sizeof(*transfer_struct_c.ts_do_debug_write_),transfer_struct_c.ts_do_debug_write__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_write_fourier_",sizeof(*transfer_struct_c.ts_do_write_fourier_),transfer_struct_c.ts_do_write_fourier__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_write_input_",sizeof(*transfer_struct_c.ts_do_write_input_),transfer_struct_c.ts_do_write_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_write_results_",sizeof(*transfer_struct_c.ts_do_write_results_),transfer_struct_c.ts_do_write_results__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_write_scenario_",sizeof(*transfer_struct_c.ts_do_write_scenario_),transfer_struct_c.ts_do_write_scenario__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_fixed_write transfer_struct_c;

  
};

// Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_fixed_inputs_c_alloc_init(struct lidort_fixed_inputs *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_inputs_c_init_only(struct lidort_fixed_inputs *transfer_struct_c, void **fortran_type_c);
  void lidort_fixed_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_fixed_inputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_fixed_inputs {
  void* bool_;
  int bool__f_byte_size;

  void* chapman_;
  int chapman__f_byte_size;

  void* cont_;
  int cont__f_byte_size;

  void* optical_;
  int optical__f_byte_size;

  void* sunrays_;
  int sunrays__f_byte_size;

  void* userval_;
  int userval__f_byte_size;

  void* write_;
  int write__f_byte_size;

  
};

// Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Fixed_Inputs : public Spurr_Fixed_Inputs_Base {
public:
  // Allocating constructor
  Lidort_Fixed_Inputs() : Spurr_Fixed_Inputs_Base() {
    lidort_fixed_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Fixed_Inputs(void* allocated_f_type_c) : Spurr_Fixed_Inputs_Base(allocated_f_type_c) {
    lidort_fixed_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Fixed_Inputs() {
    if (owns_pointer)
      lidort_fixed_inputs_c_destroy(&fortran_type_c);
  }

  Lidort_Fixed_Boolean& f_bool() {
    return *bool_;
  }

  const Lidort_Fixed_Boolean& f_bool() const {
    return *bool_;
  }

  void f_bool(Lidort_Fixed_Boolean& bool_in) {
    void* src_ptr = bool_in.fortran_type_ptr();
    void* dst_ptr = bool_->fortran_type_ptr();
    lidort_fixed_boolean_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Fixed_Chapman& chapman() {
    return *chapman_;
  }

  const Lidort_Fixed_Chapman& chapman() const {
    return *chapman_;
  }

  void chapman(Lidort_Fixed_Chapman& chapman_in) {
    void* src_ptr = chapman_in.fortran_type_ptr();
    void* dst_ptr = chapman_->fortran_type_ptr();
    lidort_fixed_chapman_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Fixed_Control& cont() {
    return *cont_;
  }

  const Lidort_Fixed_Control& cont() const {
    return *cont_;
  }

  void cont(Lidort_Fixed_Control& cont_in) {
    void* src_ptr = cont_in.fortran_type_ptr();
    void* dst_ptr = cont_->fortran_type_ptr();
    lidort_fixed_control_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Fixed_Optical& optical() {
    return *optical_;
  }

  const Lidort_Fixed_Optical& optical() const {
    return *optical_;
  }

  void optical(Lidort_Fixed_Optical& optical_in) {
    void* src_ptr = optical_in.fortran_type_ptr();
    void* dst_ptr = optical_->fortran_type_ptr();
    lidort_fixed_optical_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Fixed_Sunrays& sunrays() {
    return *sunrays_;
  }

  const Lidort_Fixed_Sunrays& sunrays() const {
    return *sunrays_;
  }

  void sunrays(Lidort_Fixed_Sunrays& sunrays_in) {
    void* src_ptr = sunrays_in.fortran_type_ptr();
    void* dst_ptr = sunrays_->fortran_type_ptr();
    lidort_fixed_sunrays_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Fixed_Uservalues& userval() {
    return *userval_;
  }

  const Lidort_Fixed_Uservalues& userval() const {
    return *userval_;
  }

  void userval(Lidort_Fixed_Uservalues& userval_in) {
    void* src_ptr = userval_in.fortran_type_ptr();
    void* dst_ptr = userval_->fortran_type_ptr();
    lidort_fixed_uservalues_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Fixed_Write& write() {
    return *write_;
  }

  const Lidort_Fixed_Write& write() const {
    return *write_;
  }

  void write(Lidort_Fixed_Write& write_in) {
    void* src_ptr = write_in.fortran_type_ptr();
    void* dst_ptr = write_->fortran_type_ptr();
    lidort_fixed_write_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Fixed_Inputs:" << std::endl
      << "   bool: " << f_bool()  << std::endl
      << "chapman: " << chapman()  << std::endl
      << "   cont: " << cont()  << std::endl
      << "optical: " << optical()  << std::endl
      << "sunrays: " << sunrays()  << std::endl
      << "userval: " << userval()  << std::endl
      << "  write: " << write()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    bool_.reset( new Lidort_Fixed_Boolean(transfer_struct_c.bool_) );
    chapman_.reset( new Lidort_Fixed_Chapman(transfer_struct_c.chapman_) );
    cont_.reset( new Lidort_Fixed_Control(transfer_struct_c.cont_) );
    optical_.reset( new Lidort_Fixed_Optical(transfer_struct_c.optical_) );
    sunrays_.reset( new Lidort_Fixed_Sunrays(transfer_struct_c.sunrays_) );
    userval_.reset( new Lidort_Fixed_Uservalues(transfer_struct_c.userval_) );
    write_.reset( new Lidort_Fixed_Write(transfer_struct_c.write_) );
    
  }

  struct lidort_fixed_inputs transfer_struct_c;

  boost::shared_ptr<Lidort_Fixed_Boolean> bool_;
  boost::shared_ptr<Lidort_Fixed_Chapman> chapman_;
  boost::shared_ptr<Lidort_Fixed_Control> cont_;
  boost::shared_ptr<Lidort_Fixed_Optical> optical_;
  boost::shared_ptr<Lidort_Fixed_Sunrays> sunrays_;
  boost::shared_ptr<Lidort_Fixed_Uservalues> userval_;
  boost::shared_ptr<Lidort_Fixed_Write> write_;
  
};

// Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_modified_boolean_c_alloc_init(struct lidort_modified_boolean *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_boolean_c_init_only(struct lidort_modified_boolean *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_boolean_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_boolean_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_boolean {
  int* ts_do_additional_mvout_;
  int ts_do_additional_mvout__f_byte_size;

  int* ts_do_all_fourier_;
  int ts_do_all_fourier__f_byte_size;

  int* ts_do_bvp_telescoping_;
  int ts_do_bvp_telescoping__f_byte_size;

  int* ts_do_chapman_function_;
  int ts_do_chapman_function__f_byte_size;

  int* ts_do_deltam_scaling_;
  int ts_do_deltam_scaling__f_byte_size;

  int* ts_do_double_convtest_;
  int ts_do_double_convtest__f_byte_size;

  int* ts_do_doublet_geometry_;
  int ts_do_doublet_geometry__f_byte_size;

  int* ts_do_external_wleave_;
  int ts_do_external_wleave__f_byte_size;

  int* ts_do_focorr_;
  int ts_do_focorr__f_byte_size;

  int* ts_do_focorr_external_;
  int ts_do_focorr_external__f_byte_size;

  int* ts_do_focorr_nadir_;
  int ts_do_focorr_nadir__f_byte_size;

  int* ts_do_focorr_outgoing_;
  int ts_do_focorr_outgoing__f_byte_size;

  int* ts_do_isotropic_only_;
  int ts_do_isotropic_only__f_byte_size;

  int* ts_do_mvout_only_;
  int ts_do_mvout_only__f_byte_size;

  int* ts_do_no_azimuth_;
  int ts_do_no_azimuth__f_byte_size;

  int* ts_do_observation_geometry_;
  int ts_do_observation_geometry__f_byte_size;

  int* ts_do_rayleigh_only_;
  int ts_do_rayleigh_only__f_byte_size;

  int* ts_do_refractive_geometry_;
  int ts_do_refractive_geometry__f_byte_size;

  int* ts_do_solar_sources_;
  int ts_do_solar_sources__f_byte_size;

  int* ts_do_solution_saving_;
  int ts_do_solution_saving__f_byte_size;

  int* ts_do_sscorr_truncation_;
  int ts_do_sscorr_truncation__f_byte_size;

  int* ts_do_sscorr_usephasfunc_;
  int ts_do_sscorr_usephasfunc__f_byte_size;

  int* ts_do_thermal_transonly_;
  int ts_do_thermal_transonly__f_byte_size;

  int* ts_do_user_streams_;
  int ts_do_user_streams__f_byte_size;

  
};

// Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Modified_Boolean : public Spurr_Modified_Boolean_Base {
public:
  // Allocating constructor
  Lidort_Modified_Boolean() : Spurr_Modified_Boolean_Base() {
    lidort_modified_boolean_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Boolean(void* allocated_f_type_c) : Spurr_Modified_Boolean_Base(allocated_f_type_c) {
    lidort_modified_boolean_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Boolean() {
    if (owns_pointer)
      lidort_modified_boolean_c_destroy(&fortran_type_c);
  }

  const bool ts_do_additional_mvout() const {
    return *transfer_struct_c.ts_do_additional_mvout_ != 0;
  }

  void ts_do_additional_mvout(const bool& ts_do_additional_mvout_in) {
    *transfer_struct_c.ts_do_additional_mvout_ = ts_do_additional_mvout_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_all_fourier() const {
    return *transfer_struct_c.ts_do_all_fourier_ != 0;
  }

  void ts_do_all_fourier(const bool& ts_do_all_fourier_in) {
    *transfer_struct_c.ts_do_all_fourier_ = ts_do_all_fourier_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_bvp_telescoping() const {
    return *transfer_struct_c.ts_do_bvp_telescoping_ != 0;
  }

  void ts_do_bvp_telescoping(const bool& ts_do_bvp_telescoping_in) {
    *transfer_struct_c.ts_do_bvp_telescoping_ = ts_do_bvp_telescoping_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_chapman_function() const {
    return *transfer_struct_c.ts_do_chapman_function_ != 0;
  }

  void ts_do_chapman_function(const bool& ts_do_chapman_function_in) {
    *transfer_struct_c.ts_do_chapman_function_ = ts_do_chapman_function_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_deltam_scaling() const {
    return *transfer_struct_c.ts_do_deltam_scaling_ != 0;
  }

  void ts_do_deltam_scaling(const bool& ts_do_deltam_scaling_in) {
    *transfer_struct_c.ts_do_deltam_scaling_ = ts_do_deltam_scaling_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_double_convtest() const {
    return *transfer_struct_c.ts_do_double_convtest_ != 0;
  }

  void ts_do_double_convtest(const bool& ts_do_double_convtest_in) {
    *transfer_struct_c.ts_do_double_convtest_ = ts_do_double_convtest_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_doublet_geometry() const {
    return *transfer_struct_c.ts_do_doublet_geometry_ != 0;
  }

  void ts_do_doublet_geometry(const bool& ts_do_doublet_geometry_in) {
    *transfer_struct_c.ts_do_doublet_geometry_ = ts_do_doublet_geometry_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_external_wleave() const {
    return *transfer_struct_c.ts_do_external_wleave_ != 0;
  }

  void ts_do_external_wleave(const bool& ts_do_external_wleave_in) {
    *transfer_struct_c.ts_do_external_wleave_ = ts_do_external_wleave_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_focorr() const {
    return *transfer_struct_c.ts_do_focorr_ != 0;
  }

  void ts_do_focorr(const bool& ts_do_focorr_in) {
    *transfer_struct_c.ts_do_focorr_ = ts_do_focorr_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_focorr_external() const {
    return *transfer_struct_c.ts_do_focorr_external_ != 0;
  }

  void ts_do_focorr_external(const bool& ts_do_focorr_external_in) {
    *transfer_struct_c.ts_do_focorr_external_ = ts_do_focorr_external_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_focorr_nadir() const {
    return *transfer_struct_c.ts_do_focorr_nadir_ != 0;
  }

  void ts_do_focorr_nadir(const bool& ts_do_focorr_nadir_in) {
    *transfer_struct_c.ts_do_focorr_nadir_ = ts_do_focorr_nadir_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_focorr_outgoing() const {
    return *transfer_struct_c.ts_do_focorr_outgoing_ != 0;
  }

  void ts_do_focorr_outgoing(const bool& ts_do_focorr_outgoing_in) {
    *transfer_struct_c.ts_do_focorr_outgoing_ = ts_do_focorr_outgoing_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_isotropic_only() const {
    return *transfer_struct_c.ts_do_isotropic_only_ != 0;
  }

  void ts_do_isotropic_only(const bool& ts_do_isotropic_only_in) {
    *transfer_struct_c.ts_do_isotropic_only_ = ts_do_isotropic_only_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_mvout_only() const {
    return *transfer_struct_c.ts_do_mvout_only_ != 0;
  }

  void ts_do_mvout_only(const bool& ts_do_mvout_only_in) {
    *transfer_struct_c.ts_do_mvout_only_ = ts_do_mvout_only_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_no_azimuth() const {
    return *transfer_struct_c.ts_do_no_azimuth_ != 0;
  }

  void ts_do_no_azimuth(const bool& ts_do_no_azimuth_in) {
    *transfer_struct_c.ts_do_no_azimuth_ = ts_do_no_azimuth_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_observation_geometry() const {
    return *transfer_struct_c.ts_do_observation_geometry_ != 0;
  }

  void ts_do_observation_geometry(const bool& ts_do_observation_geometry_in) {
    *transfer_struct_c.ts_do_observation_geometry_ = ts_do_observation_geometry_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_rayleigh_only() const {
    return *transfer_struct_c.ts_do_rayleigh_only_ != 0;
  }

  void ts_do_rayleigh_only(const bool& ts_do_rayleigh_only_in) {
    *transfer_struct_c.ts_do_rayleigh_only_ = ts_do_rayleigh_only_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_refractive_geometry() const {
    return *transfer_struct_c.ts_do_refractive_geometry_ != 0;
  }

  void ts_do_refractive_geometry(const bool& ts_do_refractive_geometry_in) {
    *transfer_struct_c.ts_do_refractive_geometry_ = ts_do_refractive_geometry_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_solar_sources() const {
    return *transfer_struct_c.ts_do_solar_sources_ != 0;
  }

  void ts_do_solar_sources(const bool& ts_do_solar_sources_in) {
    *transfer_struct_c.ts_do_solar_sources_ = ts_do_solar_sources_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_solution_saving() const {
    return *transfer_struct_c.ts_do_solution_saving_ != 0;
  }

  void ts_do_solution_saving(const bool& ts_do_solution_saving_in) {
    *transfer_struct_c.ts_do_solution_saving_ = ts_do_solution_saving_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_sscorr_truncation() const {
    return *transfer_struct_c.ts_do_sscorr_truncation_ != 0;
  }

  void ts_do_sscorr_truncation(const bool& ts_do_sscorr_truncation_in) {
    *transfer_struct_c.ts_do_sscorr_truncation_ = ts_do_sscorr_truncation_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_sscorr_usephasfunc() const {
    return *transfer_struct_c.ts_do_sscorr_usephasfunc_ != 0;
  }

  void ts_do_sscorr_usephasfunc(const bool& ts_do_sscorr_usephasfunc_in) {
    *transfer_struct_c.ts_do_sscorr_usephasfunc_ = ts_do_sscorr_usephasfunc_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_thermal_transonly() const {
    return *transfer_struct_c.ts_do_thermal_transonly_ != 0;
  }

  void ts_do_thermal_transonly(const bool& ts_do_thermal_transonly_in) {
    *transfer_struct_c.ts_do_thermal_transonly_ = ts_do_thermal_transonly_in ? FORTRAN_TRUE_INT : 0;
  }

  
  const bool ts_do_user_streams() const {
    return *transfer_struct_c.ts_do_user_streams_ != 0;
  }

  void ts_do_user_streams(const bool& ts_do_user_streams_in) {
    *transfer_struct_c.ts_do_user_streams_ = ts_do_user_streams_in ? FORTRAN_TRUE_INT : 0;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Boolean:" << std::endl
      << "    ts_do_additional_mvout: " << ts_do_additional_mvout()  << std::endl
      << "         ts_do_all_fourier: " << ts_do_all_fourier()  << std::endl
      << "     ts_do_bvp_telescoping: " << ts_do_bvp_telescoping()  << std::endl
      << "    ts_do_chapman_function: " << ts_do_chapman_function()  << std::endl
      << "      ts_do_deltam_scaling: " << ts_do_deltam_scaling()  << std::endl
      << "     ts_do_double_convtest: " << ts_do_double_convtest()  << std::endl
      << "    ts_do_doublet_geometry: " << ts_do_doublet_geometry()  << std::endl
      << "     ts_do_external_wleave: " << ts_do_external_wleave()  << std::endl
      << "              ts_do_focorr: " << ts_do_focorr()  << std::endl
      << "     ts_do_focorr_external: " << ts_do_focorr_external()  << std::endl
      << "        ts_do_focorr_nadir: " << ts_do_focorr_nadir()  << std::endl
      << "     ts_do_focorr_outgoing: " << ts_do_focorr_outgoing()  << std::endl
      << "      ts_do_isotropic_only: " << ts_do_isotropic_only()  << std::endl
      << "          ts_do_mvout_only: " << ts_do_mvout_only()  << std::endl
      << "          ts_do_no_azimuth: " << ts_do_no_azimuth()  << std::endl
      << "ts_do_observation_geometry: " << ts_do_observation_geometry()  << std::endl
      << "       ts_do_rayleigh_only: " << ts_do_rayleigh_only()  << std::endl
      << " ts_do_refractive_geometry: " << ts_do_refractive_geometry()  << std::endl
      << "       ts_do_solar_sources: " << ts_do_solar_sources()  << std::endl
      << "     ts_do_solution_saving: " << ts_do_solution_saving()  << std::endl
      << "   ts_do_sscorr_truncation: " << ts_do_sscorr_truncation()  << std::endl
      << "  ts_do_sscorr_usephasfunc: " << ts_do_sscorr_usephasfunc()  << std::endl
      << "   ts_do_thermal_transonly: " << ts_do_thermal_transonly()  << std::endl
      << "        ts_do_user_streams: " << ts_do_user_streams()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_do_additional_mvout_",sizeof(*transfer_struct_c.ts_do_additional_mvout_),transfer_struct_c.ts_do_additional_mvout__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_all_fourier_",sizeof(*transfer_struct_c.ts_do_all_fourier_),transfer_struct_c.ts_do_all_fourier__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_bvp_telescoping_",sizeof(*transfer_struct_c.ts_do_bvp_telescoping_),transfer_struct_c.ts_do_bvp_telescoping__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_chapman_function_",sizeof(*transfer_struct_c.ts_do_chapman_function_),transfer_struct_c.ts_do_chapman_function__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_deltam_scaling_",sizeof(*transfer_struct_c.ts_do_deltam_scaling_),transfer_struct_c.ts_do_deltam_scaling__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_double_convtest_",sizeof(*transfer_struct_c.ts_do_double_convtest_),transfer_struct_c.ts_do_double_convtest__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_doublet_geometry_",sizeof(*transfer_struct_c.ts_do_doublet_geometry_),transfer_struct_c.ts_do_doublet_geometry__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_external_wleave_",sizeof(*transfer_struct_c.ts_do_external_wleave_),transfer_struct_c.ts_do_external_wleave__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_focorr_",sizeof(*transfer_struct_c.ts_do_focorr_),transfer_struct_c.ts_do_focorr__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_focorr_external_",sizeof(*transfer_struct_c.ts_do_focorr_external_),transfer_struct_c.ts_do_focorr_external__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_focorr_nadir_",sizeof(*transfer_struct_c.ts_do_focorr_nadir_),transfer_struct_c.ts_do_focorr_nadir__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_focorr_outgoing_",sizeof(*transfer_struct_c.ts_do_focorr_outgoing_),transfer_struct_c.ts_do_focorr_outgoing__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_isotropic_only_",sizeof(*transfer_struct_c.ts_do_isotropic_only_),transfer_struct_c.ts_do_isotropic_only__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_mvout_only_",sizeof(*transfer_struct_c.ts_do_mvout_only_),transfer_struct_c.ts_do_mvout_only__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_no_azimuth_",sizeof(*transfer_struct_c.ts_do_no_azimuth_),transfer_struct_c.ts_do_no_azimuth__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_observation_geometry_",sizeof(*transfer_struct_c.ts_do_observation_geometry_),transfer_struct_c.ts_do_observation_geometry__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_rayleigh_only_",sizeof(*transfer_struct_c.ts_do_rayleigh_only_),transfer_struct_c.ts_do_rayleigh_only__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_refractive_geometry_",sizeof(*transfer_struct_c.ts_do_refractive_geometry_),transfer_struct_c.ts_do_refractive_geometry__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_solar_sources_",sizeof(*transfer_struct_c.ts_do_solar_sources_),transfer_struct_c.ts_do_solar_sources__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_solution_saving_",sizeof(*transfer_struct_c.ts_do_solution_saving_),transfer_struct_c.ts_do_solution_saving__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_sscorr_truncation_",sizeof(*transfer_struct_c.ts_do_sscorr_truncation_),transfer_struct_c.ts_do_sscorr_truncation__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_sscorr_usephasfunc_",sizeof(*transfer_struct_c.ts_do_sscorr_usephasfunc_),transfer_struct_c.ts_do_sscorr_usephasfunc__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_thermal_transonly_",sizeof(*transfer_struct_c.ts_do_thermal_transonly_),transfer_struct_c.ts_do_thermal_transonly__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_do_user_streams_",sizeof(*transfer_struct_c.ts_do_user_streams_),transfer_struct_c.ts_do_user_streams__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_boolean transfer_struct_c;

  
};

// Links to type: "lidort_modified_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_modified_control_c_alloc_init(struct lidort_modified_control *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_control_c_init_only(struct lidort_modified_control *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_control_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_control_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_control {
  int* ts_nmoments_input_;
  int ts_nmoments_input__f_byte_size;

  
};

// Links to type: "lidort_modified_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Modified_Control : public Spurr_Modified_Control_Base {
public:
  // Allocating constructor
  Lidort_Modified_Control() : Spurr_Modified_Control_Base() {
    lidort_modified_control_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Control(void* allocated_f_type_c) : Spurr_Modified_Control_Base(allocated_f_type_c) {
    lidort_modified_control_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Control() {
    if (owns_pointer)
      lidort_modified_control_c_destroy(&fortran_type_c);
  }

  const int& ts_nmoments_input() const {
    return *transfer_struct_c.ts_nmoments_input_;
  }

  void ts_nmoments_input(const int& ts_nmoments_input_in) {
    *transfer_struct_c.ts_nmoments_input_ = ts_nmoments_input_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Control:" << std::endl
      << "ts_nmoments_input: " << ts_nmoments_input()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_nmoments_input_",sizeof(*transfer_struct_c.ts_nmoments_input_),transfer_struct_c.ts_nmoments_input__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_control transfer_struct_c;

  
};

// Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_modified_sunrays_c_alloc_init(struct lidort_modified_sunrays *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_sunrays_c_init_only(struct lidort_modified_sunrays *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_sunrays_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_sunrays_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_sunrays {
  double* ts_beam_szas_;
  int ts_beam_szas__f_shapes[1];
  int ts_beam_szas__f_byte_size;

  int* ts_nbeams_;
  int ts_nbeams__f_byte_size;

  
};

// Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Modified_Sunrays : public Spurr_Modified_Sunrays_Base {
public:
  // Allocating constructor
  Lidort_Modified_Sunrays() : Spurr_Modified_Sunrays_Base() {
    lidort_modified_sunrays_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Sunrays(void* allocated_f_type_c) : Spurr_Modified_Sunrays_Base(allocated_f_type_c) {
    lidort_modified_sunrays_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Sunrays() {
    if (owns_pointer)
      lidort_modified_sunrays_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 1>& ts_beam_szas() const {
    return ts_beam_szas_;
  }

  void ts_beam_szas(const blitz::Array<double, 1>& ts_beam_szas_in) {
    ts_beam_szas_ = ts_beam_szas_in;
  }

  
  const int& ts_nbeams() const {
    return *transfer_struct_c.ts_nbeams_;
  }

  void ts_nbeams(const int& ts_nbeams_in) {
    *transfer_struct_c.ts_nbeams_ = ts_nbeams_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Sunrays:" << std::endl
      << "ts_beam_szas: " << std::endl << ts_beam_szas()  << std::endl
      << "   ts_nbeams: " << ts_nbeams()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_beam_szas_",sizeof(*transfer_struct_c.ts_beam_szas_),transfer_struct_c.ts_beam_szas__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_nbeams_",sizeof(*transfer_struct_c.ts_nbeams_),transfer_struct_c.ts_nbeams__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_beam_szas_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_beam_szas_,
      blitz::shape(transfer_struct_c.ts_beam_szas__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_sunrays transfer_struct_c;

  blitz::Array<double, 1> ts_beam_szas_;
  
};

// Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_modified_uservalues_c_alloc_init(struct lidort_modified_uservalues *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_uservalues_c_init_only(struct lidort_modified_uservalues *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_uservalues_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_uservalues_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_uservalues {
  double* ts_geometry_specheight_;
  int ts_geometry_specheight__f_byte_size;

  int* ts_n_user_doublets_;
  int ts_n_user_doublets__f_byte_size;

  int* ts_n_user_obsgeoms_;
  int ts_n_user_obsgeoms__f_byte_size;

  int* ts_n_user_relazms_;
  int ts_n_user_relazms__f_byte_size;

  int* ts_n_user_streams_;
  int ts_n_user_streams__f_byte_size;

  double* ts_user_angles_input_;
  int ts_user_angles_input__f_shapes[1];
  int ts_user_angles_input__f_byte_size;

  double* ts_user_doublets_;
  int ts_user_doublets__f_shapes[2];
  int ts_user_doublets__f_byte_size;

  double* ts_user_levels_;
  int ts_user_levels__f_shapes[1];
  int ts_user_levels__f_byte_size;

  double* ts_user_obsgeoms_input_;
  int ts_user_obsgeoms_input__f_shapes[2];
  int ts_user_obsgeoms_input__f_byte_size;

  double* ts_user_relazms_;
  int ts_user_relazms__f_shapes[1];
  int ts_user_relazms__f_byte_size;

  
};

// Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Modified_Uservalues : public Spurr_Modified_Uservalues_Base {
public:
  // Allocating constructor
  Lidort_Modified_Uservalues() : Spurr_Modified_Uservalues_Base() {
    lidort_modified_uservalues_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Uservalues(void* allocated_f_type_c) : Spurr_Modified_Uservalues_Base(allocated_f_type_c) {
    lidort_modified_uservalues_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Uservalues() {
    if (owns_pointer)
      lidort_modified_uservalues_c_destroy(&fortran_type_c);
  }

  const double& ts_geometry_specheight() const {
    return *transfer_struct_c.ts_geometry_specheight_;
  }

  void ts_geometry_specheight(const double& ts_geometry_specheight_in) {
    *transfer_struct_c.ts_geometry_specheight_ = ts_geometry_specheight_in;
  }

  
  const int& ts_n_user_doublets() const {
    return *transfer_struct_c.ts_n_user_doublets_;
  }

  void ts_n_user_doublets(const int& ts_n_user_doublets_in) {
    *transfer_struct_c.ts_n_user_doublets_ = ts_n_user_doublets_in;
  }

  
  const int& ts_n_user_obsgeoms() const {
    return *transfer_struct_c.ts_n_user_obsgeoms_;
  }

  void ts_n_user_obsgeoms(const int& ts_n_user_obsgeoms_in) {
    *transfer_struct_c.ts_n_user_obsgeoms_ = ts_n_user_obsgeoms_in;
  }

  
  const int& ts_n_user_relazms() const {
    return *transfer_struct_c.ts_n_user_relazms_;
  }

  void ts_n_user_relazms(const int& ts_n_user_relazms_in) {
    *transfer_struct_c.ts_n_user_relazms_ = ts_n_user_relazms_in;
  }

  
  const int& ts_n_user_streams() const {
    return *transfer_struct_c.ts_n_user_streams_;
  }

  void ts_n_user_streams(const int& ts_n_user_streams_in) {
    *transfer_struct_c.ts_n_user_streams_ = ts_n_user_streams_in;
  }

  
  const blitz::Array<double, 1>& ts_user_angles_input() const {
    return ts_user_angles_input_;
  }

  void ts_user_angles_input(const blitz::Array<double, 1>& ts_user_angles_input_in) {
    ts_user_angles_input_ = ts_user_angles_input_in;
  }

  
  const blitz::Array<double, 2>& ts_user_doublets() const {
    return ts_user_doublets_;
  }

  void ts_user_doublets(const blitz::Array<double, 2>& ts_user_doublets_in) {
    ts_user_doublets_ = ts_user_doublets_in;
  }

  
  const blitz::Array<double, 1>& ts_user_levels() const {
    return ts_user_levels_;
  }

  void ts_user_levels(const blitz::Array<double, 1>& ts_user_levels_in) {
    ts_user_levels_ = ts_user_levels_in;
  }

  
  const blitz::Array<double, 2>& ts_user_obsgeoms_input() const {
    return ts_user_obsgeoms_input_;
  }

  void ts_user_obsgeoms_input(const blitz::Array<double, 2>& ts_user_obsgeoms_input_in) {
    ts_user_obsgeoms_input_ = ts_user_obsgeoms_input_in;
  }

  
  const blitz::Array<double, 1>& ts_user_relazms() const {
    return ts_user_relazms_;
  }

  void ts_user_relazms(const blitz::Array<double, 1>& ts_user_relazms_in) {
    ts_user_relazms_ = ts_user_relazms_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Uservalues:" << std::endl
      << "ts_geometry_specheight: " << ts_geometry_specheight()  << std::endl
      << "    ts_n_user_doublets: " << ts_n_user_doublets()  << std::endl
      << "    ts_n_user_obsgeoms: " << ts_n_user_obsgeoms()  << std::endl
      << "     ts_n_user_relazms: " << ts_n_user_relazms()  << std::endl
      << "     ts_n_user_streams: " << ts_n_user_streams()  << std::endl
      << "  ts_user_angles_input: " << std::endl << ts_user_angles_input()  << std::endl
      << "      ts_user_doublets: " << std::endl << ts_user_doublets()  << std::endl
      << "        ts_user_levels: " << std::endl << ts_user_levels()  << std::endl
      << "ts_user_obsgeoms_input: " << std::endl << ts_user_obsgeoms_input()  << std::endl
      << "       ts_user_relazms: " << std::endl << ts_user_relazms()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_geometry_specheight_",sizeof(*transfer_struct_c.ts_geometry_specheight_),transfer_struct_c.ts_geometry_specheight__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_user_doublets_",sizeof(*transfer_struct_c.ts_n_user_doublets_),transfer_struct_c.ts_n_user_doublets__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_user_obsgeoms_",sizeof(*transfer_struct_c.ts_n_user_obsgeoms_),transfer_struct_c.ts_n_user_obsgeoms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_user_relazms_",sizeof(*transfer_struct_c.ts_n_user_relazms_),transfer_struct_c.ts_n_user_relazms__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_n_user_streams_",sizeof(*transfer_struct_c.ts_n_user_streams_),transfer_struct_c.ts_n_user_streams__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_angles_input_",sizeof(*transfer_struct_c.ts_user_angles_input_),transfer_struct_c.ts_user_angles_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_doublets_",sizeof(*transfer_struct_c.ts_user_doublets_),transfer_struct_c.ts_user_doublets__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_levels_",sizeof(*transfer_struct_c.ts_user_levels_),transfer_struct_c.ts_user_levels__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_obsgeoms_input_",sizeof(*transfer_struct_c.ts_user_obsgeoms_input_),transfer_struct_c.ts_user_obsgeoms_input__f_byte_size);
    BYTE_SIZE_ERROR_CHECK("ts_user_relazms_",sizeof(*transfer_struct_c.ts_user_relazms_),transfer_struct_c.ts_user_relazms__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_user_angles_input_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_user_angles_input_,
      blitz::shape(transfer_struct_c.ts_user_angles_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_user_doublets_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_user_doublets_,
      blitz::shape(transfer_struct_c.ts_user_doublets__f_shapes[0],
                   transfer_struct_c.ts_user_doublets__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_user_levels_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_user_levels_,
      blitz::shape(transfer_struct_c.ts_user_levels__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    ts_user_obsgeoms_input_.reference(blitz::Array<double, 2>(transfer_struct_c.ts_user_obsgeoms_input_,
      blitz::shape(transfer_struct_c.ts_user_obsgeoms_input__f_shapes[0],
                   transfer_struct_c.ts_user_obsgeoms_input__f_shapes[1]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<2>()));
    ts_user_relazms_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_user_relazms_,
      blitz::shape(transfer_struct_c.ts_user_relazms__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_uservalues transfer_struct_c;

  blitz::Array<double, 1> ts_user_angles_input_;
  blitz::Array<double, 2> ts_user_doublets_;
  blitz::Array<double, 1> ts_user_levels_;
  blitz::Array<double, 2> ts_user_obsgeoms_input_;
  blitz::Array<double, 1> ts_user_relazms_;
  
};

// Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_modified_chapman_c_alloc_init(struct lidort_modified_chapman *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_chapman_c_init_only(struct lidort_modified_chapman *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_chapman_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_chapman_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_chapman {
  double* ts_earth_radius_;
  int ts_earth_radius__f_byte_size;

  
};

// Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Modified_Chapman : public Spurr_Modified_Chapman_Base {
public:
  // Allocating constructor
  Lidort_Modified_Chapman() : Spurr_Modified_Chapman_Base() {
    lidort_modified_chapman_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Chapman(void* allocated_f_type_c) : Spurr_Modified_Chapman_Base(allocated_f_type_c) {
    lidort_modified_chapman_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Chapman() {
    if (owns_pointer)
      lidort_modified_chapman_c_destroy(&fortran_type_c);
  }

  const double& ts_earth_radius() const {
    return *transfer_struct_c.ts_earth_radius_;
  }

  void ts_earth_radius(const double& ts_earth_radius_in) {
    *transfer_struct_c.ts_earth_radius_ = ts_earth_radius_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Chapman:" << std::endl
      << "ts_earth_radius: " << ts_earth_radius()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_earth_radius_",sizeof(*transfer_struct_c.ts_earth_radius_),transfer_struct_c.ts_earth_radius__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_chapman transfer_struct_c;

  
};

// Links to type: "lidort_modified_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_modified_optical_c_alloc_init(struct lidort_modified_optical *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_optical_c_init_only(struct lidort_modified_optical *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_optical_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_optical_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_optical {
  double* ts_omega_total_input_;
  int ts_omega_total_input__f_shapes[1];
  int ts_omega_total_input__f_byte_size;

  
};

// Links to type: "lidort_modified_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Modified_Optical : public Spurr_Modified_Optical_Base {
public:
  // Allocating constructor
  Lidort_Modified_Optical() : Spurr_Modified_Optical_Base() {
    lidort_modified_optical_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Optical(void* allocated_f_type_c) : Spurr_Modified_Optical_Base(allocated_f_type_c) {
    lidort_modified_optical_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Optical() {
    if (owns_pointer)
      lidort_modified_optical_c_destroy(&fortran_type_c);
  }

  const blitz::Array<double, 1>& ts_omega_total_input() const {
    return ts_omega_total_input_;
  }

  void ts_omega_total_input(const blitz::Array<double, 1>& ts_omega_total_input_in) {
    ts_omega_total_input_ = ts_omega_total_input_in;
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Optical:" << std::endl
      << "ts_omega_total_input: " << std::endl << ts_omega_total_input()  << std::endl;

  }

  void check_byte_sizes() {
    BYTE_SIZE_ERROR_CHECK("ts_omega_total_input_",sizeof(*transfer_struct_c.ts_omega_total_input_),transfer_struct_c.ts_omega_total_input__f_byte_size);
    
  }

private:
  void link_blitz_arrays() {
    ts_omega_total_input_.reference(blitz::Array<double, 1>(transfer_struct_c.ts_omega_total_input_,
      blitz::shape(transfer_struct_c.ts_omega_total_input__f_shapes[0]),
      blitz::neverDeleteData, blitz::ColumnMajorArray<1>()));
    
  }

  void link_nested_types() {
    
  }

  struct lidort_modified_optical transfer_struct_c;

  blitz::Array<double, 1> ts_omega_total_input_;
  
};

// Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
extern "C" {
  void lidort_modified_inputs_c_alloc_init(struct lidort_modified_inputs *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_inputs_c_init_only(struct lidort_modified_inputs *transfer_struct_c, void **fortran_type_c);
  void lidort_modified_inputs_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void lidort_modified_inputs_c_destroy(void **fortran_type_c);
  
}

struct lidort_modified_inputs {
  void* mbool_;
  int mbool__f_byte_size;

  void* mchapman_;
  int mchapman__f_byte_size;

  void* mcont_;
  int mcont__f_byte_size;

  void* moptical_;
  int moptical__f_byte_size;

  void* msunrays_;
  int msunrays__f_byte_size;

  void* muserval_;
  int muserval__f_byte_size;

  
};

// Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.f90"
class Lidort_Modified_Inputs : public Spurr_Modified_Inputs_Base {
public:
  // Allocating constructor
  Lidort_Modified_Inputs() : Spurr_Modified_Inputs_Base() {
    lidort_modified_inputs_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  Lidort_Modified_Inputs(void* allocated_f_type_c) : Spurr_Modified_Inputs_Base(allocated_f_type_c) {
    lidort_modified_inputs_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~Lidort_Modified_Inputs() {
    if (owns_pointer)
      lidort_modified_inputs_c_destroy(&fortran_type_c);
  }

  Lidort_Modified_Boolean& mbool() {
    return *mbool_;
  }

  const Lidort_Modified_Boolean& mbool() const {
    return *mbool_;
  }

  void mbool(Lidort_Modified_Boolean& mbool_in) {
    void* src_ptr = mbool_in.fortran_type_ptr();
    void* dst_ptr = mbool_->fortran_type_ptr();
    lidort_modified_boolean_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Modified_Chapman& mchapman() {
    return *mchapman_;
  }

  const Lidort_Modified_Chapman& mchapman() const {
    return *mchapman_;
  }

  void mchapman(Lidort_Modified_Chapman& mchapman_in) {
    void* src_ptr = mchapman_in.fortran_type_ptr();
    void* dst_ptr = mchapman_->fortran_type_ptr();
    lidort_modified_chapman_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Modified_Control& mcont() {
    return *mcont_;
  }

  const Lidort_Modified_Control& mcont() const {
    return *mcont_;
  }

  void mcont(Lidort_Modified_Control& mcont_in) {
    void* src_ptr = mcont_in.fortran_type_ptr();
    void* dst_ptr = mcont_->fortran_type_ptr();
    lidort_modified_control_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Modified_Optical& moptical() {
    return *moptical_;
  }

  const Lidort_Modified_Optical& moptical() const {
    return *moptical_;
  }

  void moptical(Lidort_Modified_Optical& moptical_in) {
    void* src_ptr = moptical_in.fortran_type_ptr();
    void* dst_ptr = moptical_->fortran_type_ptr();
    lidort_modified_optical_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Modified_Sunrays& msunrays() {
    return *msunrays_;
  }

  const Lidort_Modified_Sunrays& msunrays() const {
    return *msunrays_;
  }

  void msunrays(Lidort_Modified_Sunrays& msunrays_in) {
    void* src_ptr = msunrays_in.fortran_type_ptr();
    void* dst_ptr = msunrays_->fortran_type_ptr();
    lidort_modified_sunrays_c_copy(&src_ptr, &dst_ptr);
  }

  
  Lidort_Modified_Uservalues& muserval() {
    return *muserval_;
  }

  const Lidort_Modified_Uservalues& muserval() const {
    return *muserval_;
  }

  void muserval(Lidort_Modified_Uservalues& muserval_in) {
    void* src_ptr = muserval_in.fortran_type_ptr();
    void* dst_ptr = muserval_->fortran_type_ptr();
    lidort_modified_uservalues_c_copy(&src_ptr, &dst_ptr);
  }

  
  

  
  virtual void print(std::ostream &output_stream) const {
    output_stream << "Lidort_Modified_Inputs:" << std::endl
      << "   mbool: " << mbool()  << std::endl
      << "mchapman: " << mchapman()  << std::endl
      << "   mcont: " << mcont()  << std::endl
      << "moptical: " << moptical()  << std::endl
      << "msunrays: " << msunrays()  << std::endl
      << "muserval: " << muserval()  << std::endl;

  }

  void check_byte_sizes() {
    
  }

private:
  void link_blitz_arrays() {
    
  }

  void link_nested_types() {
    mbool_.reset( new Lidort_Modified_Boolean(transfer_struct_c.mbool_) );
    mchapman_.reset( new Lidort_Modified_Chapman(transfer_struct_c.mchapman_) );
    mcont_.reset( new Lidort_Modified_Control(transfer_struct_c.mcont_) );
    moptical_.reset( new Lidort_Modified_Optical(transfer_struct_c.moptical_) );
    msunrays_.reset( new Lidort_Modified_Sunrays(transfer_struct_c.msunrays_) );
    muserval_.reset( new Lidort_Modified_Uservalues(transfer_struct_c.muserval_) );
    
  }

  struct lidort_modified_inputs transfer_struct_c;

  boost::shared_ptr<Lidort_Modified_Boolean> mbool_;
  boost::shared_ptr<Lidort_Modified_Chapman> mchapman_;
  boost::shared_ptr<Lidort_Modified_Control> mcont_;
  boost::shared_ptr<Lidort_Modified_Optical> moptical_;
  boost::shared_ptr<Lidort_Modified_Sunrays> msunrays_;
  boost::shared_ptr<Lidort_Modified_Uservalues> muserval_;
  
};



}
#endif