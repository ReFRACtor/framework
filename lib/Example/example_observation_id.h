#ifndef EXAMPLE_OBSERVATION_ID_H
#define EXAMPLE_OBSERVATION_ID_H

#include <boost/shared_ptr.hpp>
#include "hdf_file.h"
#include "observation_id.h"

namespace FullPhysics {
/****************************************************************//**
  Defines a simple mechanism by where a specific dataset from a
  HDF file is read along with some sort of identifer to obtain
  an index that should be used for extracting data.
*******************************************************************/
template<class T> class ExampleObservationId: public ObservationId {
public:
    ExampleObservationId(const boost::shared_ptr<HdfFile>& input, const std::string& dataset_name, const T& identifier);
    int data_index() const {return data_index_;}
private:
    int data_index_;
};

}
#endif
