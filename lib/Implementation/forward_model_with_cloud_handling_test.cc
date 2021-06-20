#include "forward_model_with_cloud_handling.h"
#include "pressure_with_cloud_handling.h"
#include "ground_with_cloud_handling.h"
#include "standard_forward_model.h"
#include "atmosphere_standard.h"
#include "unit_test_support.h"
#include "fp_serialize_support.h"
#include "lidort_rt.h"

using namespace FullPhysics;
using namespace blitz;

class PrintSpectrum: public Observer<boost::shared_ptr<NamedSpectrum> >
{
public:
  virtual ~PrintSpectrum() {}
  virtual void notify_update(const boost::shared_ptr<NamedSpectrum>& S)
  {
    std::cerr << S->name() << ": " << S->spectral_range().data() << "\n";
  };
  
};

BOOST_FIXTURE_TEST_SUITE(forward_model_with_cloud_handling, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // Note this currently contains a hardcoded path to the absco
  // data. Not clear how to handle this in general. We'll move to the
  // cross section version instead of absco for this particular case,
  // but it would be good to be able to handle this in some general
  // sort of way. A filename class of some kind that handles
  // environment variables? So other method?
  auto underlying_fm = serialize_read<StandardForwardModel>(test_data_dir() + "in/forward_model/omi_underlying_forward_model.xml");
  // We don't save these values, should perhaps change the serialized
  // data
  auto rt = boost::dynamic_pointer_cast<LidortRt>(underlying_fm->radiative_transfer());
  auto lid_interface = rt->rt_driver()->lidort_interface();
  lid_interface->lidort_modin().mbool().ts_do_sscorr_nadir(false);
  lid_interface->lidort_modin().mbool().ts_do_sscorr_outgoing(false);
  lid_interface->lidort_modin().mbool().ts_do_rayleigh_only(true);
  lid_interface->lidort_modin().mbool().ts_do_double_convtest(false);
  lid_interface->lidort_modin().mbool().ts_do_deltam_scaling(false);
  lid_interface->lidort_modin().mchapman().ts_earth_radius(6371.0);
  lid_interface->lidort_fixin().cont().ts_nstreams(2);
  lid_interface->lidort_modin().mcont().ts_nmoments_input(2);
  
  std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >
    cloud_handling_vector;
  auto atm = boost::dynamic_pointer_cast<AtmosphereStandard>(rt->atmosphere());
  cloud_handling_vector.push_back(boost::dynamic_pointer_cast<PressureWithCloudHandling>(atm->pressure_ptr()));
  cloud_handling_vector.push_back(boost::dynamic_pointer_cast<GroundWithCloudHandling>(atm->ground()));
  ForwardModelWithCloudHandling fm(underlying_fm,
	   boost::make_shared<CloudFractionFromState>(0.35332658886909485),
	   cloud_handling_vector);
  PrintSpectrum pspec;
  fm.add_observer(pspec);
  Spectrum rcfrac = fm.radiance(0, true);
}

BOOST_AUTO_TEST_CASE(serialization)
{
  if(!have_serialize_supported())
    return;
}

BOOST_AUTO_TEST_SUITE_END()

