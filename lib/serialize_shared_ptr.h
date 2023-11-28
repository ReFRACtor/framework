#include <iostream>
#include "swig_type_mapper_base.h"
#include <boost/serialization/shared_ptr_helper.hpp>
#ifdef SWIG_HAVE_BOOST_SERIALIZATION
#include <boost/archive/polymorphic_xml_iarchive.hpp>
#include <boost/archive/polymorphic_xml_oarchive.hpp>
#include <boost/archive/polymorphic_binary_iarchive.hpp>
#include <boost/archive/polymorphic_binary_oarchive.hpp>
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

  boost::serialization::shared_ptr_helper<boost::shared_ptr> & h =
    ar.template get_helper<shared_ptr_helper<boost::shared_ptr> >(
            shared_ptr_helper_id);
  h.reset(t,r);
  
  // This rest of this was added to our copy to support swig python directors
  boost::shared_ptr<SWIG_MAPPER_NAMESPACE::GenericObject> v2 =
    boost::dynamic_pointer_cast<SWIG_MAPPER_NAMESPACE::GenericObject>(t);
  if(v2) {
    if(SWIG_MAPPER_NAMESPACE::SwigTypeMapperBase::is_python_director(v2)) {
      SWIG_MAPPER_NAMESPACE::SwigTypeMapperBase::swig_python_director_setup(v2);
    }
  }
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
