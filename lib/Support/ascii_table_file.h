#ifndef ASCII_TABLE_FILE_H
#define ASCII_TABLE_FILE_H

#include <blitz/array.h>
#include <iostream>

#include "printable.h"

namespace FullPhysics {

/****************************************************************//**
 Reads an ASCII file that contains a simple table of values with
 each row having the same number of columns. Intended mostly for
 unit testing.
*******************************************************************/

class AsciiTableFile: public Printable<AsciiTableFile> {
public:

    //-----------------------------------------------------------------------
    /// Reads an ASCII file with simple rows and columns of doubles and
    /// save into the data attribute of the object
    //-----------------------------------------------------------------------

    AsciiTableFile(const std::string& filename, int skip_rows = 0);
    
    ~AsciiTableFile() = default;

    virtual void print(std::ostream& Os) const;

    virtual const blitz::Array<double, 2> data() const { return data_; }

protected:
    blitz::Array<double, 2> data_;
    
    virtual void read_data(std::istream& stream);
};
}

#endif
