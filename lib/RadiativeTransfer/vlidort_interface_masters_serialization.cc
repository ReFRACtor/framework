#include "vlidort_interface_masters.h"
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
// Serialization functions of VBrdf_Linsup_Masters
template<class Archive>
void VBrdf_Linsup_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VBrdf_Linsup_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VBrdf_Linsup_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void VBrdf_Linsup_Masters::load(Archive & ar, const unsigned int UNUSED(version))
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
// Serialization functions of VBrdf_Sup_Masters
template<class Archive>
void VBrdf_Sup_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VBrdf_Sup_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VBrdf_Sup_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void VBrdf_Sup_Masters::load(Archive & ar, const unsigned int UNUSED(version))
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
// Serialization functions of VLidort_Inputs
template<class Archive>
void VLidort_Inputs::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Inputs);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Inputs::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void VLidort_Inputs::load(Archive & ar, const unsigned int UNUSED(version))
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
// Serialization functions of VLidort_Masters
template<class Archive>
void VLidort_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void VLidort_Masters::load(Archive & ar, const unsigned int UNUSED(version))
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
// Serialization functions of VLidort_L_Inputs
template<class Archive>
void VLidort_L_Inputs::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_L_Inputs);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_L_Inputs::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void VLidort_L_Inputs::load(Archive & ar, const unsigned int UNUSED(version))
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
// Serialization functions of VLidort_Lcs_Masters
template<class Archive>
void VLidort_Lcs_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Lcs_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Lcs_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void VLidort_Lcs_Masters::load(Archive & ar, const unsigned int UNUSED(version))
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
// Serialization functions of VLidort_Lps_Masters
template<class Archive>
void VLidort_Lps_Masters::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Lps_Masters);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Lps_Masters::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void VLidort_Lps_Masters::load(Archive & ar, const unsigned int UNUSED(version))
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
// Serialization functions of VLidort_Vbrdf_Sup_Accessories
template<class Archive>
void VLidort_Vbrdf_Sup_Accessories::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(VLidort_Vbrdf_Sup_Accessories);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void VLidort_Vbrdf_Sup_Accessories::save(Archive & ar, const unsigned int UNUSED(version)) const
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
void VLidort_Vbrdf_Sup_Accessories::load(Archive & ar, const unsigned int UNUSED(version))
{
  TemporaryFile tf;
  std::string data;
  ar & FP_NVP(data);
  std::ofstream of(tf.file_name.string(), std::ios::binary);
  of.write(data.c_str(), data.size());
  of.close();
  read_fortran_file(tf.file_name.string());
}

FP_IMPLEMENT(VBrdf_Linsup_Masters);
FP_IMPLEMENT(VBrdf_Sup_Masters);
FP_IMPLEMENT(VLidort_Inputs);
FP_IMPLEMENT(VLidort_Masters);
FP_IMPLEMENT(VLidort_L_Inputs);
FP_IMPLEMENT(VLidort_Lcs_Masters);
FP_IMPLEMENT(VLidort_Lps_Masters);
FP_IMPLEMENT(VLidort_Vbrdf_Sup_Accessories);

#endif