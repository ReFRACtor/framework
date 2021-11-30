#include "lidort_interface_masters.h"
#include "fp_serialize_support.h"
#include <fstream>
#include <boost/filesystem.hpp>

// This was written by hand. Might be nice to include with automatic
// code generation.

using namespace FullPhysics;

// Temporary files. We can move this to a more central place if
// needed, but for now we just use this here.

namespace FullPhysics {
class TemporaryFile{
public:
  TemporaryFile() {
    file_name = boost::filesystem::unique_path();
  }
  ~TemporaryFile() {
    boost::filesystem::remove(file_name);
  }
  boost::filesystem::path file_name;
};
}

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Brdf_Lin_Sup_Masters::serialize(Archive& ar,
					     const unsigned int version)
{
  FP_GENERIC_BASE(Brdf_Lin_Sup_Masters);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Brdf_Lin_Sup_Masters::save(Archive & ar,
		    const unsigned int UNUSED(version)) const
{
  TemporaryFile tf;
  write_fortran_file(tf.file_name.string());
  std::string data;
  std::ifstream fin(tf.file_name.string(), std::ios::binary);
  auto sz = boost::filesystem::file_size(tf.file_name);
  data.resize(sz, '\0');
  fin.read(&data[0], sz);
  ar & FP_NVP(data);
}
template<class Archive>
void Brdf_Lin_Sup_Masters::load(Archive & ar,
			    const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

template<class Archive>
void Lidort_Lps_Masters::serialize(Archive& ar,
					     const unsigned int version)
{
  FP_GENERIC_BASE(Lidort_Lps_Masters);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Lidort_Lps_Masters::save(Archive & ar,
		    const unsigned int UNUSED(version)) const
{
  TemporaryFile tf;
  write_fortran_file(tf.file_name.string());
  std::string data;
  std::ifstream fin(tf.file_name.string(), std::ios::binary);
  auto sz = boost::filesystem::file_size(tf.file_name);
  data.resize(sz, '\0');
  fin.read(&data[0], sz);
  ar & FP_NVP(data);
}
template<class Archive>
void Lidort_Lps_Masters::load(Archive & ar,
			    const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

FP_IMPLEMENT(Brdf_Lin_Sup_Masters);
FP_IMPLEMENT(Lidort_Lps_Masters);
#endif
