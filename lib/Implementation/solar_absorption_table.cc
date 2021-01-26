#include "solar_absorption_table.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "linear_algebra.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SolarAbsorptionTable::serialize(Archive & ar,
			const unsigned int version)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SolarAbsorptionSpectrum)
    // table is pretty large, so we don't save it. Instead, we reread
    // the file
    // & FP_NVP(domain_unit) & FP_NVP(table)
    & FP_NVP(hdf_file_name) & FP_NVP(hdf_group);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void SolarAbsorptionTable::save(Archive & UNUSED(ar),
				const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void SolarAbsorptionTable::load(Archive & UNUSED(ar),
				const unsigned int UNUSED(version)) 
{
  HdfFile f(hdf_file_name);
  init(f, hdf_group);
}

FP_IMPLEMENT(SolarAbsorptionTable);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SolarAbsorptionTable, SolarAbsorptionSpectrum)
.def(luabind::constructor<const HdfFile&,
     const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read the given file for the solar absorption spectrum.
//-----------------------------------------------------------------------

SolarAbsorptionTable::SolarAbsorptionTable(
const HdfFile& F,
const std::string& Hdf_group)
{
  init(F, Hdf_group);
}

void SolarAbsorptionTable::init(const HdfFile& F, const std::string& Hdf_group)
{
  hdf_file_name = F.file_name();
  hdf_group = Hdf_group;
  ArrayWithUnit<double, 1> sdom = F.read_field_with_unit<double,1>
    (Hdf_group + "/wavenumber");
  ArrayWithUnit<double, 1> srange = F.read_field_with_unit<double,1>
    (Hdf_group + "/spectrum");
  if(!srange.units.is_commensurate(Unit("dimensionless")))
    throw Exception("Solar absorption spectrum units need to be dimensionless");
  if(sdom.value.rows() != srange.value.rows())
    throw Exception("wavenumber and spectrum need to be the same size");
  domain_unit = sdom.units;
  Array<double, 1> x = to_c_order(sdom.value);
  Array<double, 1> y = to_c_order(srange.value);
  table = LinearInterpolate<double, double>
    (x.dataFirst(), x.dataFirst() + x.rows(), y.dataFirst());
}

// See base class for description.
Spectrum SolarAbsorptionTable::solar_absorption_spectrum(
const SpectralDomain& spec_domain) const
{
  Array<double, 1> wv = spec_domain.convert_wave(domain_unit);

  Array<double, 1> res;
  res.resize(wv.shape());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = table(wv(i));
  return Spectrum(spec_domain, SpectralRange(res, Unit("dimensionless")));
}

void SolarAbsorptionTable::print(std::ostream& Os) const
{ 
  Os << "SolarAbsorptionTable\n"
     << "  Hdf file name: " << hdf_file_name << "\n"
     << "  Hdf group:     " << hdf_group << "\n";
}
