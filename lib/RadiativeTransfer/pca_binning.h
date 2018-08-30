#ifndef PCA_BINNING_H
#define PCA_BINNING_H

#include "atmosphere_oco.h"
#include "spectral_domain.h"

namespace FullPhysics {

/****************************************************************//**
 *******************************************************************/

class PCABinning {
public:
    PCABinning(const boost::shared_ptr<AtmosphereOco>& atm) : atmosphere(atm) {};

    void compute_bins(const SpectralDomain& spec_domain, int chanel_index, std::string primary_absorber);

    const blitz::Array<int, 1> ncnt() const { return ncnt_; }
    const blitz::Array<int, 1> index() const { return index_; }
    const blitz::Array<int, 1> bin() const { return bin_; }

private:
    boost::shared_ptr<AtmosphereOco> atmosphere;

    blitz::Array<int, 1> ncnt_;
    blitz::Array<int, 1> index_;
    blitz::Array<int, 1> bin_;
 
};

}

#endif
