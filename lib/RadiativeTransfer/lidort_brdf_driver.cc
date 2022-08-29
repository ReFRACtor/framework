#include "lidort_brdf_driver.h"

#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "ostream_pad.h"
#include "fe_disable_exception.h"
#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION

template<class Archive>
void LidortBrdfDriver::serialize(Archive & ar,
                                 const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpurrBrdfDriver)
    & FP_NVP_(nstream) & FP_NVP_(nmoment)
    & FP_NVP_(brdf_interface);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void LidortBrdfDriver::save(Archive & UNUSED(a),
                    const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void LidortBrdfDriver::load(Archive & UNUSED(ar),
                            const unsigned int UNUSED(version))
{
  Brdf_Sup_Inputs& brdf_inputs = brdf_interface_->brdf_sup_in();
  brdf_params.reference( brdf_inputs.bs_brdf_parameters() );
  brdf_factors.reference( brdf_inputs.bs_brdf_factors() );
}

FP_IMPLEMENT(LidortBrdfDriver);

#endif

//=======================================================================
// LidortBrdfInterface
//=======================================================================

//-----------------------------------------------------------------------
/// Initialize Lidort BRDF interface
//-----------------------------------------------------------------------

LidortBrdfDriver::LidortBrdfDriver(int nstream, int nmoment)
  : nstream_(nstream),
    nmoment_(nmoment)
{
  init();
}

void LidortBrdfDriver::init()
{
  brdf_interface_.reset( new Brdf_Lin_Sup_Masters() );

  Brdf_Sup_Inputs& brdf_inputs = brdf_interface_->brdf_sup_in();

  // Only use 1 beam meaning only one set of sza, azm
  brdf_inputs.bs_nbeams(1);
  brdf_inputs.bs_n_user_streams(1);
  brdf_inputs.bs_n_user_relazms(1);

  // This MUST be consistent with streams used for
  // LIDORT RT calculation
  brdf_inputs.bs_nstreams(nstream_);

  // Recommended value from LIDORT manual
  // Number of quadtrature streams for BRDF calculation
  brdf_inputs.bs_nstreams_brdf(50);

  brdf_params.reference( brdf_inputs.bs_brdf_parameters() );
  brdf_factors.reference( brdf_inputs.bs_brdf_factors() );
}

void LidortBrdfDriver::setup_geometry(double sza, double azm, double zen)
{

  Brdf_Sup_Inputs& brdf_inputs = brdf_interface_->brdf_sup_in();

  // Solar zenith angles (degrees) [0,90]
  Array<double, 1> brdf_sza( brdf_inputs.bs_beam_szas() );
  brdf_sza(0) = sza;

  // User-defined relative angles (in degrees) for
  // off-quadrature output.
  Array<double, 1> brdf_azm( brdf_inputs.bs_user_relazms() );
  brdf_azm(0) = azm;

  // User-defined viewing zenith angles (in degrees) for
  // off quadrature output.
  Array<double, 1> brdf_zen( brdf_inputs.bs_user_angles_input() );
  brdf_zen(0) = zen;
}

void LidortBrdfDriver::calculate_brdf() const
{
  // Lidort may cause floating point exceptions when doing setup. This
  // is because it may copy garbage value, which are never used. By
  // chance the garbage values may cause a overflow. We
  // suspend floating point exceptions when doing setup
  FeDisableException disable_fp;

  // Process BRDF inputs
  bool do_debug_restoration = false;
  brdf_interface_->run(do_debug_restoration, nmoment_);
}

int LidortBrdfDriver::n_brdf_kernels() const
{
  return brdf_interface_->brdf_sup_in().bs_n_brdf_kernels();
}

void LidortBrdfDriver::n_brdf_kernels(const int n_kernels)
{
  brdf_interface_->brdf_sup_in().bs_n_brdf_kernels(n_kernels);
}

int LidortBrdfDriver::n_kernel_factor_wfs() const {
  return brdf_interface_->brdf_linsup_in().bs_n_kernel_factor_wfs();
}

void LidortBrdfDriver::n_kernel_factor_wfs(const int n_factors) {
  brdf_interface_->brdf_linsup_in().bs_n_kernel_factor_wfs(n_factors);
}

int LidortBrdfDriver::n_kernel_params_wfs() const {
  return brdf_interface_->brdf_linsup_in().bs_n_kernel_params_wfs();
}

void LidortBrdfDriver::n_kernel_params_wfs(const int n_params) {
  brdf_interface_->brdf_linsup_in().bs_n_kernel_params_wfs(n_params);
}

int LidortBrdfDriver::n_surface_wfs() const {
  return brdf_interface_->brdf_linsup_in().bs_n_surface_wfs();
}

void LidortBrdfDriver::n_surface_wfs(const int n_wfs) {
  brdf_interface_->brdf_linsup_in().bs_n_surface_wfs(n_wfs);
}

bool LidortBrdfDriver::do_kparams_derivs(const int kernel_index) const
{
  return brdf_interface_->brdf_linsup_in().bs_do_kparams_derivs()(kernel_index);
}

void LidortBrdfDriver::do_kparams_derivs(const int kernel_index, const bool do_kparams)
{
  Brdf_Linsup_Inputs& brdf_lin_inputs = brdf_interface_->brdf_linsup_in();
  Array<bool, 1> do_kparams_derivs( brdf_lin_inputs.bs_do_kparams_derivs() );
  do_kparams_derivs(kernel_index) = do_kparams;
  brdf_lin_inputs.bs_do_kparams_derivs(do_kparams_derivs);
}

bool LidortBrdfDriver::do_shadow_effect() const {
  return brdf_interface_->brdf_sup_in().bs_do_shadow_effect();
}

void LidortBrdfDriver::do_shadow_effect(const bool do_shadow) const {
  brdf_interface_->brdf_sup_in().bs_do_shadow_effect(do_shadow);
}

void LidortBrdfDriver::initialize_kernel_parameters(const int kernel_index,
                                                    const int which_brdf,
                                                    const bool lambertian_flag,
                                                    const int n_brdf_parameters,
                                                    const bool do_factor_wfs,
                                                    const blitz::Array<bool, 1>& do_params_wfs)
{
  Brdf_Sup_Inputs& brdf_inputs = brdf_interface_->brdf_sup_in();
  Brdf_Linsup_Inputs& brdf_lin_inputs = brdf_interface_->brdf_linsup_in();

  Array<int, 1> bs_which_brdf( brdf_inputs.bs_which_brdf() );
  bs_which_brdf(kernel_index) = which_brdf;

  Array<bool, 1> bs_lambertian_flag = brdf_inputs.bs_lambertian_kernel_flag();
  bs_lambertian_flag(kernel_index) = lambertian_flag;
  brdf_inputs.bs_lambertian_kernel_flag(bs_lambertian_flag);

  Array<int, 1> bs_n_brdf_parameters( brdf_inputs.bs_n_brdf_parameters() );
  bs_n_brdf_parameters(kernel_index) = n_brdf_parameters;

  Array<bool, 1> bs_do_factor_wfs( brdf_lin_inputs.bs_do_kernel_factor_wfs() );
  bs_do_factor_wfs(kernel_index) = do_factor_wfs;
  brdf_lin_inputs.bs_do_kernel_factor_wfs(bs_do_factor_wfs);

  Array<bool, 2> bs_do_params_wfs( brdf_lin_inputs.bs_do_kernel_params_wfs() );
  bs_do_params_wfs(kernel_index, Range(0, do_params_wfs.rows()-1)) = do_params_wfs;
  brdf_lin_inputs.bs_do_kernel_params_wfs(bs_do_params_wfs);
}

