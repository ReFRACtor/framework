#include "error_analysis.h"
#include "absorber_absco.h"
#include <new> // std::bad_alloc
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(ErrorAnalysis)
.def(luabind::constructor<const boost::shared_ptr<ConnorSolver>&,
     const boost::shared_ptr<RtAtmosphere>&,
     const boost::shared_ptr<ForwardModel>&,
     const boost::shared_ptr<Observation>&>())
.def(luabind::constructor<const boost::shared_ptr<MaxAPosteriori>&,
     const boost::shared_ptr<RtAtmosphere>&,
     const boost::shared_ptr<ForwardModel>&,
     const boost::shared_ptr<Observation>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// This is variation of ErrorAnalysis that takes a general
/// RtAtmosphere. This will *fail* for anything other than a
/// AtmosphereStandard, however it is convenient to have this for use with
/// Lua (which has much more limited knowledge of the class structure).
//-----------------------------------------------------------------------

ErrorAnalysis::ErrorAnalysis(const boost::shared_ptr<ConnorSolver>& Solver,
                             const boost::shared_ptr<RtAtmosphere>& Atm,
                             const boost::shared_ptr<ForwardModel>& Fm,
                             const boost::shared_ptr<Observation>& inst_meas)
    : solver(Solver), atm(boost::dynamic_pointer_cast<AtmosphereStandard>(Atm)),
      fm(Fm), meas(inst_meas)
{
}

//-----------------------------------------------------------------------
/// This is variation of ErrorAnalysis that takes a general
/// RtAtmosphere. This will *fail* for anything other than a
/// AtmosphereStandard, however it is convenient to have this for use with
/// Lua (which has much more limited knowledge of the class structure).
//-----------------------------------------------------------------------

ErrorAnalysis::ErrorAnalysis(const boost::shared_ptr<MaxAPosteriori>& Max_a_posteriori,
                             const boost::shared_ptr<RtAtmosphere>& Atm,
                             const boost::shared_ptr<ForwardModel>& Fm,
                             const boost::shared_ptr<Observation>& inst_meas)
    : max_a_posteriori(Max_a_posteriori),
      atm(boost::dynamic_pointer_cast<AtmosphereStandard>(Atm)),
      fm(Fm), meas(inst_meas)
{
}

//-----------------------------------------------------------------------
/// The number of spectral bands associated with forward model.
//-----------------------------------------------------------------------

int ErrorAnalysis::num_channels() const
{
    return fm->num_channels();
}

//-----------------------------------------------------------------------
/// The residual is the difference between modeled and measured radiances
//-----------------------------------------------------------------------

blitz::Array<double, 1> ErrorAnalysis::residual() const
{
    if(solver) {
        return solver->residual();
    }
    else {
        return max_a_posteriori->model_measure_diff();
    }
}

//-----------------------------------------------------------------------
/// Modeled radiance.
//-----------------------------------------------------------------------

blitz::Array<double, 1> ErrorAnalysis::modeled_radiance() const
{
    if(residual().rows() == 0) {
        return blitz::Array<double, 1>(0);
    }

    return blitz::Array<double, 1>
           (residual() + meas->radiance_all().spectral_range().data());
}

//-----------------------------------------------------------------------
/// Calculate an approximation to the size of the continuum signal
/// where there is no significant atmosphere absorption. We
/// approximate this by finding the 10 highest radiance values and
/// averaging them.
//-----------------------------------------------------------------------

double ErrorAnalysis::signal_level(int band) const
{
    FeDisableException disable_fp;
    const int nrad = 10;
    SpectralRange rad(meas->radiance(band).spectral_range());

    if(rad.data().rows() == 0) {
        return 0;
    }

    Array<double, 1> r(rad.data().copy());
    std::sort(r.data(), r.data() + r.rows()); // Min to max value
    r.reverseSelf(firstDim);       // Now max to min value
    Range r2(0, std::min(nrad - 1, r.rows() - 1));
    return sum(r(r2) / r2.length());
}

//-----------------------------------------------------------------------
/// Helper class for sort done in noise.
//-----------------------------------------------------------------------
// Don't have Doxygen document this class.
/// @cond
class DataRow {
public:
    DataRow(int R, double V) : row(R), value(V) {}
    int row;
    double value;
};

class DataRowCompare {
public:
    bool operator()(const DataRow& D1, const DataRow& D2) const
    {
        return D1.value > D2.value;
    }
};

/// @endcond

//-----------------------------------------------------------------------
/// Calculate an approximation to the size noise of the continuum signal
/// where there is no significant atmosphere absorption. We
/// approximate this by finding the 10 highest radiance values and
/// averaging their noise.
//-----------------------------------------------------------------------

double ErrorAnalysis::noise_level(int band) const
{
    FeDisableException disable_fp;
    const int nrad = 10;
    SpectralRange rad(meas->radiance(band).spectral_range());

    if(rad.data().rows() < nrad || rad.uncertainty().rows() < nrad) {
        return 0;
    }

    std::vector<DataRow> dr;
    dr.reserve(rad.data().rows());

    for(int i = 0; i < rad.data().rows(); ++i) {
        dr.push_back(DataRow(i, rad.data()(i)));
    }

    std::sort(dr.begin(), dr.end(), DataRowCompare()); // Max to min value
    double sum_noise = 0;

    for(int i = 0; i < nrad; ++i) {
        sum_noise += rad.uncertainty()(dr[i].row);
    }

    return sum_noise / nrad;
}

//-----------------------------------------------------------------------
///
//-----------------------------------------------------------------------

blitz::Array<double, 2> ErrorAnalysis::aposteriori_covariance() const
{
    if(solver) {
        return solver->aposteriori_covariance();
    }
    else {
        return max_a_posteriori->a_posteriori_covariance();
    }
}

//-----------------------------------------------------------------------
///
//-----------------------------------------------------------------------

blitz::Array<double, 2> ErrorAnalysis::apriori_covariance() const
{
    if(solver) {
        return solver->apriori_covariance();
    }
    else {
        return max_a_posteriori->a_priori_cov();
    }
}

//-----------------------------------------------------------------------
///
//-----------------------------------------------------------------------

blitz::Array<double, 2> ErrorAnalysis::averaging_kernel() const
{
    if(solver) {
        return solver->averaging_kernel();
    }
    else {
        return max_a_posteriori->averaging_kernel();
    }
}


//-----------------------------------------------------------------------
/// Return the sum of the squares of the residual for the given band.
//-----------------------------------------------------------------------

double ErrorAnalysis::residual_sum_sq(int Band) const
{
    FeDisableException disable_fp;
    boost::optional<blitz::Range> pr = fm->stacked_pixel_range(Band);

    if(!pr) {
        return 0;
    }

    blitz::Array<double, 1> res(residual());

    if(res.rows() == 0) {
        return 0;
    }

    return sum(res(*pr) * res(*pr));
}

//-----------------------------------------------------------------------
/// Return the residual mean square for the O2 band.
//-----------------------------------------------------------------------

double ErrorAnalysis::residual_mean_sq(int Band) const
{
    FeDisableException disable_fp;
    boost::optional<blitz::Range> pr = fm->stacked_pixel_range(Band);

    if(!pr) {
        return 0;
    }

    return sqrt(residual_sum_sq(Band) / pr->length());
}

//-----------------------------------------------------------------------
/// Return the reduced chisq for band
//-----------------------------------------------------------------------

double ErrorAnalysis::reduced_chisq(int Band) const
{
    FeDisableException disable_fp;
    boost::optional<blitz::Range> pr = fm->stacked_pixel_range(Band);

    if(!pr) {
        return 0;
    }

    if(solver) {
        if(solver->residual().rows() == 0) {
            return 0;
        }

        return chisq_measure_norm(solver->residual()(*pr),
                                  solver->residual_covariance_diagonal()(*pr));
    }
    else {
        blitz::Array<double, 1> res(max_a_posteriori->uncert_weighted_model_measure_diff());

        if(!res.rows()) {
            return 0;
        }

        return sum(res(*pr) * res(*pr)) / pr->length();
    }
}

//-----------------------------------------------------------------------
/// Return the relative residual mean square for the given band.
//-----------------------------------------------------------------------

double ErrorAnalysis::relative_residual_mean_sq(int Band) const
{
    FeDisableException disable_fp;
    boost::optional<blitz::Range> pr = fm->stacked_pixel_range(Band);

    if(!pr) {
        return 0;
    }

    double result = residual_mean_sq(Band);
    return result ? (result / signal_level(Band)) : 0;
}

//-----------------------------------------------------------------------
/// return chisq_measure_norm for the given data.
//-----------------------------------------------------------------------

double ErrorAnalysis::chisq_measure_norm(const blitz::Array<double, 1>& Residual,
                                         const blitz::Array<double, 1>& Residual_cov_diag) const
{
    FeDisableException disable_fp;

    if (residual().rows() == 0) {
        return 0;
    }

    return sum(Residual * Residual / Residual_cov_diag) / Residual.rows();
}

//-----------------------------------------------------------------------
/// Calculate the degrees of freedom for the full state vector. This
/// is just the trace of the averaging kernel.
///
/// \todo ATB reference?
//-----------------------------------------------------------------------

double ErrorAnalysis::degrees_of_freedom() const
{
    FeDisableException disable_fp;
    return sum(averaging_kernel()(i1, i1));
}
