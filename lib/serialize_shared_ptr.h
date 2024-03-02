#include <iostream>
#include "swig_type_mapper_base.h"
#include <boost/serialization/shared_ptr_helper.hpp>
#ifdef SWIG_HAVE_BOOST_SERIALIZATION
#include <boost/archive/polymorphic_iarchive.hpp>
#include <boost/archive/polymorphic_oarchive.hpp>
#endif

//-----------------------------------------------------------------------
/// This is largely a duplicate of boost/serialization/shared_ptr.hpp.
/// We need to have our own copy to add support for handling SWIG
/// directors.
//-----------------------------------------------------------------------

namespace boost {
  namespace serialization{
template<class T>
struct version< ::boost::shared_ptr< T > > {
  typedef mpl::integral_c_tag tag;
  typedef mpl::int_<1> type;
  static const int value = type::value;
};
// don't track shared pointers
template<class T>
struct tracking_level< ::boost::shared_ptr< T > > {
  typedef mpl::integral_c_tag tag;
  typedef mpl::int_< ::boost::serialization::track_never> type;
  static const int value = type::value;
};

void * const shared_ptr_helper_id = 0;

struct null_deleter {
    void operator()(void const *) const {}
};

//-----------------------------------------------------------------------
/// MMS    
/// We actually want to modify shared_ptr_helper, but this is a bit
/// involved code and copying it seemed like it would be fragile. So
/// instead we add a special version of the shared_ptr, but just adds
/// its own version of reset and other functions used by
/// shared_ptr_helper. This pretty much is just a normal shared_ptr,
/// but adds special handling for the pointer passed to reset.
//-----------------------------------------------------------------------

template<class T> class our_shared_ptr : public boost::shared_ptr<T>
{
public:
  typedef typename boost::shared_ptr< T >::element_type element_type;
  template< class Y >
  our_shared_ptr(our_shared_ptr<Y> const & r, element_type * p )
    : boost::shared_ptr<T>(r, p)
  {
  }
  template< class Y >
  our_shared_ptr(our_shared_ptr<Y> const & r)
    : boost::shared_ptr<T>(r)
  {
  }
  our_shared_ptr() {}
  void reset() BOOST_SP_NOEXCEPT
  {
    boost::shared_ptr<T>::reset();
  }
  void reset(T* r)
  {
    // We'll put in a proper deleter shortly.
    boost::shared_ptr<T> t2(r, null_deleter());
    boost::shared_ptr<SWIG_MAPPER_NAMESPACE::GenericObject> v2 =
      boost::dynamic_pointer_cast<SWIG_MAPPER_NAMESPACE::GenericObject>(t2);
    if(v2 && SWIG_MAPPER_NAMESPACE::SwigTypeMapperBase::is_python_director(v2)) {
      boost::shared_ptr<T> t = boost::dynamic_pointer_cast<T>(SWIG_MAPPER_NAMESPACE::SwigTypeMapperBase::swig_python_director_setup(v2));
      boost::shared_ptr<T>::operator=(t);
    } else {
      boost::shared_ptr<T>::reset(r);
    }      
  }
};
  
template<class Archive, class T>
inline void save(
		 Archive & ar,
		 const boost::shared_ptr< T > &t,
		 const unsigned int /* file_version */
		 ){
  // The most common cause of trapping here would be serializing
  // something like shared_ptr<int>.  This occurs because int
  // is never tracked by default.  Wrap int in a trackable type
  BOOST_STATIC_ASSERT((tracking_level< T >::value != track_never));
  const T * t_ptr = t.get();
  ar << boost::serialization::make_nvp("px", t_ptr);
}

template<class Archive, class T>
inline void load(
    Archive & ar,
    boost::shared_ptr< T > &t,
    const unsigned int /*file_version*/
){
  // The most common cause of trapping here would be serializing
  // something like shared_ptr<int>.  This occurs because int
  // is never tracked by default.  Wrap int in a trackable type
  BOOST_STATIC_ASSERT((tracking_level< T >::value != track_never));
  T* r;
  ar >> boost::serialization::make_nvp("px", r);
  /// MMS, use our_shared_ptr to get the proper handling of
  /// PythonRefPtrCleanup in place.
  boost::serialization::shared_ptr_helper<our_shared_ptr> & h =
    ar.template get_helper<shared_ptr_helper<our_shared_ptr> >(
            shared_ptr_helper_id);
  our_shared_ptr<T> t2;
  h.reset(t2,r);
  t = t2;
}
    
template<class Archive, class T>
inline void serialize(
    Archive & ar,
    boost::shared_ptr< T > &t,
    const unsigned int file_version
){
    // correct shared_ptr serialization depends upon object tracking
    // being used.
    BOOST_STATIC_ASSERT(
        boost::serialization::tracking_level< T >::value
        != boost::serialization::track_never
    );
    boost::serialization::split_free(ar, t, file_version);
}
    
  }
}
