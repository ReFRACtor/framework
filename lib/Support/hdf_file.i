// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "hdf_file.h"
%}
%base_import(generic_object)
%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::HdfFile);

namespace FullPhysics {
class HdfFile : public GenericObject {
public:
    enum Mode {READ, CREATE, READ_WRITE};
    
    HdfFile(const std::string& Fname, Mode M = READ);
    void close();
    
    bool has_object(const std::string& Objname) const;
    bool has_attribute(const std::string& Aname) const;
  
    %python_attribute(file_name, std::string)
    %python_attribute(mode, Mode)
    static bool is_hdf(const std::string& Fname);
  
    std::string print_to_string();
  
    %extend {
      void write_double_1d(const std::string& fname, const blitz::Array<double, 1>& D)
          { $self->write_field(fname, D); }
      void write_double_2d(const std::string& fname, const blitz::Array<double, 2>& D)
          { $self->write_field(fname, D); }
      void write_double_3d(const std::string& fname, const blitz::Array<double, 3>& D)
          { $self->write_field(fname, D); }
      void write_double_4d(const std::string& fname, const blitz::Array<double, 4>& D)
          { $self->write_field(fname, D); }
      blitz::Array<double, 1> read_double_1d(const std::string& fname)
          { return $self->read_field<double, 1>(fname); }
      blitz::Array<double, 2> read_double_2d(const std::string& fname)
          { return $self->read_field<double, 2>(fname); }
      blitz::Array<double, 3> read_double_3d(const std::string& fname)
          { return $self->read_field<double, 3>(fname); }
      blitz::Array<double, 4> read_double_4d(const std::string& fname)
          { return $self->read_field<double, 4>(fname); }
      blitz::Array<int, 1> read_int_1d(const std::string& fname)
          { return $self->read_field<int, 1>(fname); }
      blitz::Array<int, 2> read_int_2d(const std::string& fname)
          { return $self->read_field<int, 2>(fname); }
      blitz::Array<int, 3> read_int_3d(const std::string& fname)
          { return $self->read_field<int, 3>(fname); }
      blitz::Array<int, 4> read_int_4d(const std::string& fname)
          { return $self->read_field<int, 4>(fname); }
      ArrayWithUnit<double, 1> read_double_with_unit_1d(const std::string& fname) 
          { return $self->read_field_with_unit<double, 1>(fname); }
      ArrayWithUnit<double, 2> read_double_with_unit_2d(const std::string& fname) 
          { return $self->read_field_with_unit<double, 2>(fname); }
      ArrayWithUnit<double, 3> read_double_with_unit_3d(const std::string& fname) 
          { return $self->read_field_with_unit<double, 3>(fname); }
      ArrayWithUnit<double, 4> read_double_with_unit_4d(const std::string& fname) 
          { return $self->read_field_with_unit<double, 4>(fname); }

      std::string read_attribute_string(const std::string& aname)
          { return $self->read_attribute<std::string>(aname); } 
  
      template<class T, int D> ArrayWithUnit<T, D> read_field_with_unit(const std::string& Dataname) const;
    }

  %pickle_serialization();
};

}
