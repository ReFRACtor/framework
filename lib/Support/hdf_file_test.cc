#include "unit_test_support.h"
#include "hdf_file.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(hdf_file, GlobalFixture)

BOOST_AUTO_TEST_CASE(read)
{
  BOOST_CHECK(HdfFile::is_hdf(test_data_dir() + "in/common/l1b_example_data.h5"));
  BOOST_CHECK(!HdfFile::is_hdf(test_data_dir() + "in/solver/connor_converged"));
  HdfFile h(test_data_dir() + "in/common/l1b_example_data.h5");
  BOOST_CHECK(h.has_object("/Level1b/stokes_coefficient"));
  BOOST_CHECK(!h.has_object("/Level1b/fake_field"));
  Array<double, 2> st_expect(3,4);
  st_expect = 
    0.5, 0.499993, -0.00268922, 0,
    0.5, 0.499993, -0.00269793, 0,
    0.5, 0.499993, -0.0025686, 0;
  Array<double, 3> st = h.read_field<double, 3>("/Level1b/stokes_coefficient");
  Array<double, 2> stsub = st(0, Range::all(), Range::all());
  BOOST_CHECK_MATRIX_CLOSE_TOL(stsub, st_expect, 1e-5);
  TinyVector<int, 3> start, size;
  start = 0, 0, 0;
  size = 1, st.cols(), st.extent(thirdDim);
  Array<double, 3> st2 = h.read_field<double, 3>("/Level1b/stokes_coefficient", start, size);
  BOOST_CHECK_EQUAL(st2.rows(), 1);
  BOOST_CHECK_EQUAL(st2.cols(), 3);
  BOOST_CHECK_EQUAL(st2.extent(thirdDim), 4);
  Array<double, 2> stsub2 = st2(0, Range::all(), Range::all());
  BOOST_CHECK_MATRIX_CLOSE_TOL(stsub2, st_expect, 1e-5);

  BOOST_CHECK_EQUAL(h.read_attribute<std::string>("Level1b/altitude/Units"), "Meters");

  ArrayWithUnit<double, 2> alt = h.read_field_with_unit<double, 2>("/Level1b/altitude");
  BOOST_CHECK(alt.units.is_commensurate(units::m));
}

BOOST_AUTO_TEST_CASE(write)
{
  HdfFile h("hdf_file_write.h5", HdfFile::CREATE);
  add_file_to_cleanup("hdf_file_write.h5");
  Array<double, 2> d(2,3);
  d = 1, 2, 3,
    4, 5, 6;
  Array<int, 3> d2(2,3,4);
  d2 = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24;

  Array<std::string, 2> st_var(2, 3);
  st_var = "s1", "s2", "s3",
    "t1", "t2", "t3";

  // Gather c pointers from strings above
  Array<const char*, 2> st_fixed(st_var.shape());
  for (int i = 0; i < st_fixed.rows(); i++) {
    for (int j = 0; j < st_fixed.cols(); j++) {
      st_fixed(i,j) = st_var(i,j).c_str();
    }
  }

  h.write_field("/TestGroup/TestSubGroup/Test", d);
  h.write_field("/TestGroup/TestSubGroup2/Test2", d2);
  h.write_field("/TestGroup/TestSubGroup2/TestVarString", st_var);
  h.write_field("/TestGroup/TestSubGroup2/TestFixedString", st_fixed);
  h.write_field("/TestGroup/TestSubGroup2/Scalar", 5.2);
  h.write_attribute("/TestGroup/TestSubGroup/Test/att1", 4);
  h.write_attribute("/TestGroup/TestSubGroup/att2", d);
  h.write_attribute("/TestGroup2/TestSubGroup3/att3", 4);
  h.write_attribute("/TestGroup4/TestSubGroup/att4", "hi there");
  h.dimension_metadata("Dim1", "This is a description", 2);
  h.dimension_metadata("Dim2", "This is a description", 3);
  h.shape_metadata("Arr1", "Dim1");
  h.shape_metadata("Arr2", "Dim1", "Dim2");
  std::vector<std::string> svec;
  svec.push_back(std::string("string1"));
  svec.push_back(std::string("string2"));
  h.write_attribute("/TestGroup4/att5", svec);
  h.close();
  HdfFile hread("hdf_file_write.h5");
  Array<double, 2> dread = 
    hread.read_field<double, 2>("/TestGroup/TestSubGroup/Test");
  Array<int, 3> d2read = 
    hread.read_field<int, 3>("/TestGroup/TestSubGroup2/Test2");
  Array<std::string, 2> stread_var = 
    hread.read_field<std::string, 2>("/TestGroup/TestSubGroup2/TestVarString");
  Array<std::string, 2> stread_fix = 
    hread.read_field<std::string, 2>("/TestGroup/TestSubGroup2/TestFixedString");
  double sread = hread.read_field<double>("/TestGroup/TestSubGroup2/Scalar");
  int att1read = hread.read_attribute<int>("/TestGroup/TestSubGroup/Test/att1");
  Array<double, 2> att2read =
    hread.read_attribute<double, 2>("/TestGroup/TestSubGroup/att2");
  BOOST_CHECK_MATRIX_CLOSE(dread, d);
  BOOST_CHECK_MATRIX_CLOSE(d2read,  d2);
  BOOST_CHECK_EQUAL(stread_var(0,0), "s1");
  BOOST_CHECK_EQUAL(stread_var(0,1), "s2");
  BOOST_CHECK_EQUAL(stread_var(0,2), "s3");
  BOOST_CHECK_EQUAL(stread_var(1,0), "t1");
  BOOST_CHECK_EQUAL(stread_var(1,1), "t2");
  BOOST_CHECK_EQUAL(stread_var(1,2), "t3");
  BOOST_CHECK_EQUAL(stread_fix(0,0), "s1");
  BOOST_CHECK_EQUAL(stread_fix(0,1), "s2");
  BOOST_CHECK_EQUAL(stread_fix(0,2), "s3");
  BOOST_CHECK_EQUAL(stread_fix(1,0), "t1");
  BOOST_CHECK_EQUAL(stread_fix(1,1), "t2");
  BOOST_CHECK_EQUAL(stread_fix(1,2), "t3");
  BOOST_CHECK_CLOSE(sread, 5.2, 1e-5);
  BOOST_CHECK_EQUAL(att1read, 4);
  BOOST_CHECK_MATRIX_CLOSE(att2read, d);
  BOOST_CHECK_EQUAL(hread.read_attribute<int>("/TestGroup2/TestSubGroup3/att3"),
                    4);
  BOOST_CHECK_EQUAL(hread.read_attribute<std::string>
                    ("/TestGroup4/TestSubGroup/att4"), "hi there");
  std::vector<std::string> att5read = 
    hread.read_attribute<std::vector<std::string> >("/TestGroup4/att5");
  BOOST_CHECK_EQUAL((int) att5read.size(), 2);
  BOOST_CHECK_EQUAL(att5read[0], "string1");
  BOOST_CHECK_EQUAL(att5read[1], "string2");
}

BOOST_AUTO_TEST_CASE(bad_data)
{
  BOOST_CHECK_THROW(HdfFile h(test_data_dir() + "bad_file"), 
                    FullPhysics::Exception);
}

BOOST_AUTO_TEST_SUITE_END()
