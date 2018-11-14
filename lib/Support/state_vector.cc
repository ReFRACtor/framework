#include "state_vector.h"
#include "spectrum_effect.h"
#include "fp_exception.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
#include "rt_atmosphere.h"
#include "instrument.h"
#include "instrument_correction.h"
#include "stokes_coefficient.h"
#include "absorber_vmr.h"
#include "aerosol_optical.h"
#include "aerosol_extinction.h"
#include "aerosol_property.h"
#include "sample_grid.h"

void state_vector_add_observer_instrument(StateVector& Sv, Instrument& inst)
{
    Sv.add_observer(inst);
}

void state_vector_add_observer_atm(StateVector& Sv, RtAtmosphere& atm)
{
    Sv.add_observer(atm);
}

void state_vector_add_observer_scoeff(StateVector& Sv, StokesCoefficient& Scoef)
{
    Sv.add_observer(Scoef);
}

void state_vector_add_observer_spec_effect(StateVector& Sv, std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& se)
{
    BOOST_FOREACH(std::vector<boost::shared_ptr<SpectrumEffect> >& i, se) {
        BOOST_FOREACH(boost::shared_ptr<SpectrumEffect>& j, i)
        Sv.add_observer(*j);
    }
}

void state_vector_add_observer_inst_corr(StateVector& Sv, std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >& se)
{
    BOOST_FOREACH(std::vector<boost::shared_ptr<InstrumentCorrection> >& i, se) {
        BOOST_FOREACH(boost::shared_ptr<InstrumentCorrection>& j, i)
        Sv.add_observer(*j);
    }
}

void state_vector_add_observer_absorber_vmr(StateVector& Sv, AbsorberVmr& avmr)
{
    Sv.add_observer(avmr);
}

void state_vector_add_observer_aerosol_optical(StateVector& Sv, AerosolOptical& aopt)
{
    Sv.add_observer(aopt);
}

void state_vector_add_observer_aerosol_extinction(StateVector& Sv, AerosolExtinction& aext)
{
    Sv.add_observer(aext);
}

void state_vector_add_observer_aerosol_property(StateVector& Sv, AerosolProperty& aprop)
{
    Sv.add_observer(aprop);
}

void state_vector_add_observer_pressure(StateVector& Sv, Pressure& press)
{
    Sv.add_observer(press);
}

void state_vector_add_observer_temperature(StateVector& Sv, Temperature& temp)
{
    Sv.add_observer(temp);
}

void state_vector_add_observer_ground(StateVector& Sv, Ground& ground)
{
    Sv.add_observer(ground);
}


void state_vector_add_observer_disp(StateVector& Sv, SampleGrid& disp)
{
    Sv.add_observer(disp);
}

void state_vector_print_names(StateVector& Sv)
{
    BOOST_FOREACH(const std::string & s, Sv.state_vector_name())
    std::cout << s << "\n";
}
typedef void (StateVector::*us1)(const blitz::Array<double, 1>&);
typedef void (StateVector::*us2)(const blitz::Array<double, 1>&, const blitz::Array<double, 2>&);
REGISTER_LUA_CLASS(StateVector)
.def(luabind::constructor<>())
.def("add_observer", &state_vector_add_observer_instrument)
.def("add_observer", &state_vector_add_observer_atm)
.def("add_observer", &state_vector_add_observer_scoeff)
.def("add_observer", &state_vector_add_observer_spec_effect)
.def("add_observer", &state_vector_add_observer_inst_corr)
.def("add_observer", &state_vector_add_observer_absorber_vmr)
.def("add_observer", &state_vector_add_observer_aerosol_optical)
.def("add_observer", &state_vector_add_observer_aerosol_extinction)
.def("add_observer", &state_vector_add_observer_aerosol_property)
.def("add_observer", &state_vector_add_observer_pressure)
.def("add_observer", &state_vector_add_observer_temperature)
.def("add_observer", &state_vector_add_observer_ground)
.def("add_observer", &state_vector_add_observer_disp)
.def("update_state", ((us1) &StateVector::update_state))
.def("update_state", ((us2) &StateVector::update_state))
.def("print_names", &state_vector_print_names)
.def("state", &StateVector::state)
.def("state_covariance", &StateVector::state_covariance)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Update the state vector. This version sets the covariance to a
/// dummy identity matrix that is the same size as state().
//-----------------------------------------------------------------------

void StateVector::update_state(const blitz::Array<double, 1>& X)
{
    x_.resize(X.rows(), X.rows());
    x_.value() = X;
    x_.jacobian() = 0;

    for(int i = 0; i < x_.rows(); ++i) {
        x_.jacobian()(i, i) = 1;
    }

    // Only set dummy covariance values when it is empty or the state vector size change. Don't
    // overwrite anything that may have been written by the other overloaded update_state method
    if (cov_.rows() != x_.rows() or cov_.cols() == x_.rows()) {
        cov_.resize(x_.rows(), x_.rows());
        cov_ = 0;

        for(int i = 0; i < x_.rows(); ++i) {
            cov_(i, i) = 1;
        }
    }

    notify_update_do(*this);
}

//-----------------------------------------------------------------------
/// Update the state vector and covariance.
//-----------------------------------------------------------------------

void StateVector::update_state(const blitz::Array<double, 1>& X,
                               const blitz::Array<double, 2>& Cov)
{
    if(X.rows() != Cov.rows() ||
            X.rows() != Cov.cols()) {
        throw Exception("X and Cov need to be the same size when updating the StateVector");
    }

    x_.resize(X.rows(), X.rows());
    x_.value() = X;
    x_.jacobian() = 0;

    for(int i = 0; i < x_.rows(); ++i) {
        x_.jacobian()(i, i) = 1;
    }

    cov_.reference(Cov.copy());
    notify_update_do(*this);
}


//-----------------------------------------------------------------------
/// Return a Array of boolean values. The value (i) is true if the state
/// vector element X(i) is being used. This can be used to determine
/// parameters that are being ignored, e.g. the number of active
/// levels in an Aerosol is less that the size of the state vector for it.
//-----------------------------------------------------------------------

blitz::Array<bool, 1> StateVector::used_flag() const
{
    blitz::Array<bool, 1> res(state().rows());
    res = false;
    BOOST_FOREACH(const boost::weak_ptr<Observer<StateVector> >& t, olist) {
        boost::shared_ptr<StateVectorObserver> t2 =
            boost::dynamic_pointer_cast<StateVectorObserver>(t.lock());

        if(t2) {
            t2->mark_used(*this, res);
        }
    }
    return res;
}

//-----------------------------------------------------------------------
/// Return name of each state vector element.
//-----------------------------------------------------------------------

Array<std::string, 1> StateVector::state_vector_name() const
{
    Array<std::string, 1> res(std::max(state().rows(), observer_claimed_size()));

    for(int i = 0; i < res.rows(); ++i) {
        res(i) = "State vector " + boost::lexical_cast<std::string>(i + 1);
    }

    BOOST_FOREACH(const boost::weak_ptr<Observer<StateVector> >& t, olist) {
        boost::shared_ptr<StateVectorObserver> t2 =
            boost::dynamic_pointer_cast<StateVectorObserver>(t.lock());

        if(t2) {
            t2->state_vector_name(*this, res);
        }
    }
    return res;
}

void StateVector::print(std::ostream& Os) const
{
    const static int sv_num_width = 17;
    Array<std::string, 1> svname = state_vector_name();

    for(int i = 0; i < std::max(svname.rows(), state().rows()); ++i) {
        Os << std::setprecision(sv_num_width - 7)
           << std::setw(sv_num_width);

        if(i < state().rows()) {
            Os << state()(i);
        } else {
            Os << "Past end state vector";
        }

        Os << "  ";

        if(i < svname.rows()) {
            Os << svname(i);
        } else {
            Os << "Unlabeled row " << i;
        }

        Os << std::endl;

    }
}
