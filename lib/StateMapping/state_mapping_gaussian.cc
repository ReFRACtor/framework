#include "state_mapping_gaussian.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateMappingGaussian::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateMapping)
    & FP_NVP(map_name) & FP_NVP_(min_desired) & FP_NVP(linear_total)
    & FP_NVP(press);
}

FP_IMPLEMENT(StateMappingGaussian);
#endif

  
//-----------------------------------------------------------------------
/// Constructor.
/// \param Press Pressure object used for defining pressure levels
/// \param Linear_AOD Specifies whether the first gaussian parameter
///   that represents the desired total is in linear space if true,
///   logarithmic space if false.
//-----------------------------------------------------------------------

StateMappingGaussian::StateMappingGaussian
(const boost::shared_ptr<Pressure>& in_press,
 bool Linear_Total,
 double Min_Desired)
  : min_desired_(Min_Desired), linear_total(Linear_Total), press(in_press)
{
  if (linear_total)
    map_name = "gaussian_linear";
  else
    map_name = "gaussian_log";
}

//-----------------------------------------------------------------------
/// Total aerosol optical depth of the extinction values in component.
//-----------------------------------------------------------------------

AutoDerivative<double> StateMappingGaussian::total_optical_depth
(ArrayAd<double, 1> component) const
{
  ArrayAd<double, 1> pressure_grid = press->pressure_grid().value;
  AutoDerivative<double> tot_component = 0.0;
  for(int layer_idx = 0; layer_idx < press->number_layer(); layer_idx++) {
    AutoDerivative<double> delta_press =
      (pressure_grid(layer_idx + 1) - pressure_grid(layer_idx)) / 2.0;
    tot_component +=
      (delta_press * (component(layer_idx) + component(layer_idx + 1) ));
  }
  return tot_component;
}

//-----------------------------------------------------------------------
/// Calculation of forward model view of coeffs with mapping applied
//-----------------------------------------------------------------------

ArrayAd<double, 1> StateMappingGaussian::mapped_state
(const ArrayAd<double, 1>& updated_coeff) const
{
  ArrayAd<double, 1> component(press->number_level(), updated_coeff.number_variable());
  AutoDerivative<double> desired_val;
  
  if (linear_total)
    desired_val = updated_coeff(0);
  else
    desired_val = std::exp(updated_coeff(0));
  
  // Don't let component value go lower than a minimum value. Not clear if this
  // is actually what we want to do, see ticket #2252 for this
  // suggestion. We might instead want to use a constrained solver.
  if(desired_val < min_desired_)
    desired_val = min_desired_;
  
  int ngaussians = int((updated_coeff.rows() - 1) / 2);
  
  ArrayAd<double, 1> pressure_grid = press->pressure_grid().value;
  AutoDerivative<double> surface_press = press->surface_pressure().value;
  
  component.resize(press->number_level(), updated_coeff.number_variable());
  component.resize_number_variable(updated_coeff.number_variable());
  component.jacobian() = 0;
  
  for(int g_idx = 0; g_idx < ngaussians; g_idx++) {
    AutoDerivative<double> p0 = updated_coeff(g_idx * 2 + 1);
    AutoDerivative<double> sigma = updated_coeff(g_idx * 2 + 2);
    
    for(int lev = 0; lev < component.rows(); lev++) {
      AutoDerivative<double> p = pressure_grid(lev) / surface_press;
      AutoDerivative<double> g_eval =
        std::exp( -1 * (p - p0) * (p - p0) / (2 * sigma * sigma) );

      // Because its not that easy to init ArrayAd from python to 0.0
      if (g_idx == 0) {
        component(lev) = g_eval;
      } else {
        component(lev) = component(lev) + g_eval;
      }
    }
  }

  AutoDerivative<double> scaling_N;
   
  // Protect against a divide by zero
  if (total_optical_depth(component).value() != 0.0)
    scaling_N = desired_val / total_optical_depth(component);
  else
    scaling_N = 0.0;

  for(int lev = 0; lev < component.rows(); lev++)
    component(lev) = component(lev) * scaling_N;
  return component;
}
