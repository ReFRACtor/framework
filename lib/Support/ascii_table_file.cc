#include <fstream>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "ascii_table_file.h"

using namespace FullPhysics;
using namespace blitz;

AsciiTableFile::AsciiTableFile(const std::string& filename, int skip_rows)
{
    ifstream stream(filename);
    std::string line;

    // Skip any rows
    for(int skip_idx = 0; skip_idx < skip_rows; skip_rows++) {
        std::getline(stream, line);
    }

    // Read the remaining contents as data
    read_data(stream);

    stream.close();
}
    
void AsciiTableFile::read_data(std::istream& stream) {

    std::vector<Array<double, 1> > row_data;
    std::string line;

    // Read unknown number of rows
    int num_cols = 0;

    double value;
    while(std::getline(stream, line)) {
        // Remove white space
        boost::trim(line);

        // Skip rows with comments or which are completely empty
        if(line.size() == 0 || boost::find_first(line, "#")) {
            continue;
        }

        // Count number of columns in first row
        if(num_cols <= 0) {
            std::stringstream col_count_ss(line);
            while(col_count_ss >> value) num_cols++;
        }

        Array<double, 1> col_data(num_cols);

        std::stringstream line_ss(line);
        int col_idx = 0;
        while(line_ss >> value) col_data(col_idx++) = value;

        row_data.push_back(col_data);
    }

    // Copy columns into class attribute once know the number of rows
    data_.resize(row_data.size(), num_cols);

    for(int row_idx = 0; row_idx < data_.rows(); row_idx++) {
        data_(row_idx, Range::all()) = row_data[row_idx];
    }

}

void AsciiTableFile::print(std::ostream& Os) const
{
    Os << "AsciiTableFile: " << std::endl
       << "    data size: " << data_.rows() << " x " << data_.cols() << std::endl;
}
