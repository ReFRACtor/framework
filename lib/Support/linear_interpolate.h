#ifndef LINEAR_INTERPOLATE_H
#define LINEAR_INTERPOLATE_H
#include "printable.h"
#include "fp_exception.h"
#include "auto_derivative.h"
#include <map>
#include <boost/shared_ptr.hpp>

namespace FullPhysics {

template<class TX, class TY> class InterpolatePoint;
template<> class InterpolatePoint<AutoDerivative<double>, 
				  AutoDerivative<double> >;

/****************************************************************//**
  Interface for interpolating values
*******************************************************************/

template<class TX, class TY> class InterpolatePoint {
public:
  virtual const TX& x_min() const = 0;
  virtual const TX& x_max() const = 0;
  virtual TY operator()(const TX& x) const = 0;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  This just returns the same y value, always
*******************************************************************/

template<class TX, class TY> class Return1Point;
template<> class Return1Point<AutoDerivative<double>, AutoDerivative<double> >;

template<class TX, class TY> class Return1Point : 
    public InterpolatePoint<TX, TY>,
    public Printable<Return1Point<TX, TY> > {
public:
  Return1Point(const TX& x0, const TY& y0)
    : x0_(x0), y0_(y0) 
  {}
  const TX& x_min() const {return x0_; }
  const TX& x_max() const {return x0_; }
  TY operator()(const TX& UNUSED(x)) const { return y0_; }
  void print(std::ostream& Os) const { Os << "Return1Point"; }
private:
  TX x0_;
  TY y0_;
  Return1Point() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};


template<class TX, class TY> class LinearInterpolate2Point;
template<> class LinearInterpolate2Point<AutoDerivative<double>, 
					 AutoDerivative<double> >;

/****************************************************************//**
  This does linear interpolate between two points.
*******************************************************************/

template<class TX, class TY> class LinearInterpolate2Point :
    public InterpolatePoint<TX, TY>,
    public Printable<LinearInterpolate2Point<TX, TY> > {
public:
  LinearInterpolate2Point(const TX& x0, const TY& y0, const TX& x1, 
			  const TY& y1)
    : x0_(x0), x1_(x1), y0_(y0), delta_y0_((y1 - y0) / (x1 - x0))
  {
    range_max_check(x0, x1);
  }
  const TX& x_min() const {return x0_; }
  const TX& x_max() const {return x1_; }
  TY operator()(const TX& x) const { return TY(y0_ + delta_y0_ * (x - x0_)); }
  void print(std::ostream& Os) const { Os << "LinearInterpolate2Point"; }
private:
  TX x0_;
  TX x1_;
  TY y0_;
  TY delta_y0_;
  LinearInterpolate2Point() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

template<class TX, class TY> class LinearInterpolate;
template<> class LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >;
/****************************************************************//**
  This class takes a set of points and values, and linearly
  interpolates between those values. 

  This works for any type TX of x and TY of y. TX is typically a
  double or float, but it might also be something like
  AutoDerivative<double>. TX has the same constraints as a key to
  std::map, it needs to be ordered by the "<" operator.

  If you request a value for X outside of the range originally passed
  in, then the behavior depends on the value passed in the
  constructor. We can extrapolate, "clip" (i.e., use the minimum value
  for x < min, the maximum value for x > max), or trigger an
  exception. 
*******************************************************************/
template<class TX, class TY> class LinearInterpolate : 
    public Printable<LinearInterpolate<TX, TY> > {
public:
  enum BehaviorOutOfRange 
    {OUT_OF_RANGE_EXTRAPOLATE = 0, OUT_OF_RANGE_CLIP, OUT_OF_RANGE_ERROR};
  LinearInterpolate() {}

//-----------------------------------------------------------------------
/// Constructor. This takes iterators to give x and y, along with the 
/// behavior when called with out of range data.
//-----------------------------------------------------------------------

  template<class I1, class I2> LinearInterpolate(I1 xstart, I1 xend,
      I2 ystart, BehaviorOutOfRange Out_of_range = OUT_OF_RANGE_EXTRAPOLATE)
    : out_of_range(Out_of_range)
  {
    if(std::distance(xstart, xend) == 0) {
      // If no values then we can not do any interpolation, 
      // error in operator if attempted
      return;
    } else if(std::distance(xstart, xend) == 1) {
      // If only one value available then just return that one point
      // no need for interpolation
      const TX& x0 = *xstart;
      const TY& y0 = *ystart;

      inter[x0] = boost::shared_ptr<InterpolatePoint<TX, TY> >
	(new Return1Point<TX, TY>(x0, y0));

    } else while(xstart != xend) {
      // For two more points we can use linear interpolation
      const TX& x0 = *xstart++;
      const TY& y0 = *ystart++;
      if(xstart == xend)
	break;
      const TX& x1 = *xstart;
      const TY& y1 = *ystart;
      if(x1 <= x0)
	throw Exception("X needs to be sorted");
      
      inter[x1] = boost::shared_ptr<InterpolatePoint<TX, TY> >
	(new LinearInterpolate2Point<TX, TY>(x0, y0, x1, y1));
    }
  }
  TY operator()(const TX& x) const
  { 
    if(inter.empty()) {
      throw Exception("Can not interpolate since interpolation map is empty.");
    }

    typedef typename 
      std::map<TX, boost::shared_ptr<InterpolatePoint<TX, TY> > >::
      const_iterator Itype;
    Itype i = inter.lower_bound(x);
    if(i == inter.end()) { // Are we past the end?
      switch(out_of_range) {
      case OUT_OF_RANGE_ERROR:
	{
	  std::stringstream err_msg;
	  err_msg << "Linear interpolation point "
		  << x << " is past the end of interpolation range: "
		  << "[" << inter.begin()->first
		  << ", " << inter.rbegin()->first << "]";
	  throw Exception(err_msg.str());
	}
      case OUT_OF_RANGE_EXTRAPOLATE:
	return (*inter.rbegin()->second)(x);
      case OUT_OF_RANGE_CLIP:
	return (*inter.rbegin()->second)(inter.rbegin()->second->x_max());
      default:
	throw Exception("Unknown out_of_range value");
      }
    } else {
      if(x >= i->second->x_min())
	return (*i->second)(x);
      switch(out_of_range) {
      case OUT_OF_RANGE_ERROR:
	{
	  std::stringstream err_msg;
	  err_msg << "Linear interpolation point "
		  << x << " is outside of interpolation range: "
		  << "[" << inter.begin()->first
		  << ", " << inter.rbegin()->first << "]";
	  throw Exception(err_msg.str());
	}
      case OUT_OF_RANGE_EXTRAPOLATE:
	return (*i->second)(x);
      case OUT_OF_RANGE_CLIP:
	return (*i->second)(i->second->x_min());
      default:
	throw Exception("Unknown out_of_range value");
      }
    }
  }
  void print(std::ostream& Os) const { Os << "LinearInterpolate"; }
  static std::string name() { return "linear linear"; }
private:
  std::map<TX, boost::shared_ptr<InterpolatePoint<TX, TY> > > inter;
  BehaviorOutOfRange out_of_range;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

typedef InterpolatePoint<double, double>  InterpolatePoint_double_double;
typedef Return1Point<double, double>  Return1Point_double_double;
typedef LinearInterpolate2Point<double, double>  LinearInterpolate2Point_double_double;
typedef LinearInterpolate<double, double> LinearInterpolate_double_double;
typedef InterpolatePoint<double, AutoDerivative<double> >  InterpolatePoint_double_auto_derivative_double;
typedef Return1Point<double, AutoDerivative<double> >  Return1Point_double_auto_derivative_double;
typedef LinearInterpolate2Point<double, AutoDerivative<double> > LinearInterpolate2Point_double_auto_derivative_double;
typedef LinearInterpolate<double, AutoDerivative<double> > LinearInterpolate_double_auto_derivative_double;
}
// Specialization for autoderivative
#include "linear_interpolate_ad.h"

FP_EXPORT_KEY(InterpolatePoint_double_double);
FP_EXPORT_KEY(Return1Point_double_double);
FP_EXPORT_KEY(LinearInterpolate2Point_double_double);
FP_EXPORT_KEY(LinearInterpolate_double_double);

// I don't think these work correctly. The serialization is also
// pretty slow. We could probably come up with a more optimized way
// to serialize these, but this is only used in a handful of places
// and the data is easy enough to regenerate from other things we are
// already serializing. Comment out for now so we don't accidentally
// use this. We can uncomment and sort this out if it becomes an issue.
// FP_EXPORT_KEY(InterpolatePoint_double_auto_derivative_double);
// FP_EXPORT_KEY(Return1Point_double_auto_derivative_double);
// FP_EXPORT_KEY(LinearInterpolate2Point_double_auto_derivative_double);
// FP_EXPORT_KEY(LinearInterpolate_double_auto_derivative_double);
#endif

