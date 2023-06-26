#include "state_mapping_basis_matrix.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(state_mapping_basis_matrix, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  blitz::Array<double, 2> basis_matrix(3,2);
  basis_matrix =
    1., 0.,
    0.5, 0.5,
    0,1;
  StateMappingBasisMatrix m(basis_matrix);
  BOOST_CHECK_MATRIX_CLOSE_TOL(m.basis_matrix(), basis_matrix, 1e-7);
  blitz::Array<double, 2> inv_expt(2,3);
  inv_expt =
    5.0/6, 1.0/3, -1.0/6,
    -1.0/6,  1.0/3,  5.0/6;
  BOOST_CHECK_MATRIX_CLOSE_TOL(m.inverse_basis_matrix(), inv_expt, 1e-7);
  blitz::Array<double, 1> ret_grid(2);
  ret_grid = 1,2;
  blitz::Array<double, 1> map_expect(3);
  map_expect = 1,1.5,2;
  BOOST_CHECK_MATRIX_CLOSE_TOL(m.mapped_state(ret_grid).value(), map_expect,
			       1e-7);
  BOOST_CHECK_MATRIX_CLOSE_TOL(m.retrieval_state(map_expect).value(), ret_grid,
			       1e-7);
}

BOOST_AUTO_TEST_CASE(serialize)
{
  if(!have_serialize_supported())
    return;
  blitz::Array<double, 2> basis_matrix(3,2);
  basis_matrix =
    1., 0.,
    0.5, 0.5,
    0,1;
  blitz::Array<double, 2> inv_expt(2,3);
  inv_expt =
    5.0/6, 1.0/3, -1.0/6,
    -1.0/6,  1.0/3,  5.0/6;
  auto m = boost::make_shared<StateMappingBasisMatrix>(basis_matrix);
  std::string d = serialize_write_string(m);
  if(false)
    std::cerr << d;
  auto mr = serialize_read_string<StateMappingBasisMatrix>(d);
  BOOST_CHECK_MATRIX_CLOSE_TOL(mr->basis_matrix(), basis_matrix, 1e-7);
  BOOST_CHECK_MATRIX_CLOSE_TOL(mr->inverse_basis_matrix(), inv_expt, 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()
