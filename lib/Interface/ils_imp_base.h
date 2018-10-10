#ifndef ILS_IMP_BASE_H
#define ILS_IMP_BASE_H

#include "double_with_unit.h"
#include "sample_grid.h"
#include "ils.h"

namespace FullPhysics {

/****************************************************************//**
 Base class to inherit from to simplify the task of implementing 
 most of the necessary Ils behavior. Inheriting classes need to 
 define the apply_ils methods.
*******************************************************************/

class IlsImpBase : public Ils, public Observer<SampleGrid> {
public:
    //-----------------------------------------------------------------------
    /// Constructor.
    //-----------------------------------------------------------------------
    IlsImpBase(const boost::shared_ptr<SampleGrid>& Sample_grid, const DoubleWithUnit& Ils_half_width)
    : sample_grid_(Sample_grid), ils_half_width_(Ils_half_width)
    {
        sample_grid_->add_observer(*this);
    }

    virtual ~IlsImpBase() = default;

    virtual void notify_update(const SampleGrid& D)
    {
        notify_update_do(*this);
    }

    virtual blitz::Array<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const blitz::Array<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const = 0;

    virtual ArrayAd<double, 1> apply_ils
    (const blitz::Array<double, 1>& High_resolution_wave_number,
     const ArrayAd<double, 1>& High_resolution_radiance,
     const std::vector<int>& Pixel_list) const = 0;

    virtual void print(std::ostream& Os) const
    {
        Os << "IlsImpBase";
    }

    virtual const std::string band_name() const
    {
        return desc_band_name_;
    }
    virtual const std::string hdf_band_name() const
    {
        return hdf_band_name_;
    }
    virtual const SpectralDomain pixel_grid() const
    {
        return sample_grid_->pixel_grid();
    }
    virtual const DoubleWithUnit ils_half_width() const
    {
        return ils_half_width_;
    }
    virtual void ils_half_width(const DoubleWithUnit& half_width)
    {
        ils_half_width_ = half_width;
    }

    //-----------------------------------------------------------------------
    /// Underlying SampleGrid object
    //-----------------------------------------------------------------------
    boost::shared_ptr<SampleGrid> sample_grid() const {return sample_grid_; }

    virtual boost::shared_ptr<Ils> clone() const = 0;

private:
    boost::shared_ptr<SampleGrid> sample_grid_;
    DoubleWithUnit ils_half_width_;
    std::string desc_band_name_, hdf_band_name_;
};
}
#endif
