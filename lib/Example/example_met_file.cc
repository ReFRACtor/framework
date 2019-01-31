#include <boost/make_shared.hpp>

#include "example_met_file.h"
#include "example_observation_id.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ExampleMetFile, Meteorology)
.def(luabind::constructor<const boost::shared_ptr<HdfFile>&, const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// \param input_file HDF file with meteological data
/// \param observation_id Identifier string that is present in the observation_ids dataset
//-----------------------------------------------------------------------

ExampleMetFile::ExampleMetFile(const boost::shared_ptr<HdfFile>& input_file, const std::string& observation_id)
: input(input_file)
{
    ExampleObservationId<std::string> obs_id(input_file, "observation_ids", observation_id);
    data_index = obs_id.data_index();
}

ExampleMetFile::ExampleMetFile(const std::string& input_filename, const std::string& observation_id)
: input(boost::make_shared<HdfFile>(input_filename))
{
    ExampleObservationId<std::string> obs_id(input, "observation_ids", observation_id);
    data_index = obs_id.data_index();
}

//-----------------------------------------------------------------------
/// Read a field where a single number is expected to be returned
//-----------------------------------------------------------------------

double ExampleMetFile::read_scalar(const std::string& field) const
{
    return input->read_field<double, 1>(field, TinyVector<int, 1>(data_index), TinyVector<int, 1>(1))(0);
}

//-----------------------------------------------------------------------
/// Read a field and by its data index
//-----------------------------------------------------------------------

Array<double, 1> ExampleMetFile::read_array(const std::string& Field) const
{
    TinyVector<int, 2> sz = input->read_shape<2>(Field);
    Array<double, 2> raw = input->read_field<double, 2>(Field, TinyVector<int, 2>(data_index, 0), TinyVector<int, 2>(1, sz[1]));
    Array<double, 1> value(raw.extent(secondDim));

    value = raw(0, Range::all());

    return value;
}
