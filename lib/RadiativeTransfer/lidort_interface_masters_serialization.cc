#include "lidort_interface_masters.h"
#include "fp_serialize_support.h"
#include <fstream>
#include <boost/filesystem.hpp>

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

////////
// Serialization functions of Brdf_Lin_Sup_Masters
template<class Archive>
void Brdf_Lin_Sup_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Brdf_Lin_Sup_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Brdf_Lin_Sup_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void Brdf_Lin_Sup_Masters::load(Archive & ar, const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

////////
// Serialization functions of Brdf_Sup_Masters
template<class Archive>
void Brdf_Sup_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Brdf_Sup_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Brdf_Sup_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void Brdf_Sup_Masters::load(Archive & ar, const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

////////
// Serialization functions of Lidort_Inputs
template<class Archive>
void Lidort_Inputs::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Lidort_Inputs);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Lidort_Inputs::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void Lidort_Inputs::load(Archive & ar, const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

////////
// Serialization functions of Lidort_Masters
template<class Archive>
void Lidort_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Lidort_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Lidort_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void Lidort_Masters::load(Archive & ar, const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

////////
// Serialization functions of Lidort_L_Inputs
template<class Archive>
void Lidort_L_Inputs::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Lidort_L_Inputs);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Lidort_L_Inputs::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void Lidort_L_Inputs::load(Archive & ar, const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

////////
// Serialization functions of Lidort_Lcs_Masters
template<class Archive>
void Lidort_Lcs_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Lidort_Lcs_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Lidort_Lcs_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void Lidort_Lcs_Masters::load(Archive & ar, const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

////////
// Serialization functions of Lidort_Lps_Masters
template<class Archive>
void Lidort_Lps_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Lidort_Lps_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Lidort_Lps_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void Lidort_Lps_Masters::load(Archive & ar, const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

////////
// Serialization functions of Lidort_Brdf_Sup_Accessories
template<class Archive>
void Lidort_Brdf_Sup_Accessories::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Lidort_Brdf_Sup_Accessories);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Lidort_Brdf_Sup_Accessories::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void Lidort_Brdf_Sup_Accessories::load(Archive & ar, const unsigned int UNUSED(version))
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
FP_IMPLEMENT(Brdf_Sup_Masters);
FP_IMPLEMENT(Lidort_Inputs);
FP_IMPLEMENT(Lidort_Masters);
FP_IMPLEMENT(Lidort_L_Inputs);
FP_IMPLEMENT(Lidort_Lcs_Masters);
FP_IMPLEMENT(Lidort_Lps_Masters);
FP_IMPLEMENT(Lidort_Brdf_Sup_Accessories);

#endif