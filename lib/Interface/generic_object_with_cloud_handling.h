#ifndef GENERIC_OBJECT_WITH_CLOUD_HANDLING_H
#define GENERIC_OBJECT_WITH_CLOUD_HANDLING_H
#include <printable.h>

namespace FullPhysics {
//-----------------------------------------------------------------------
/// GenericObject that as a do_cloud variable. This is pretty oddly
/// specific, but we have this so ForwardModelWithCloudFraction can
/// just have a list of objects it needs to toggle rather than somehow
/// knowing all the types.  
//-----------------------------------------------------------------------

class GenericObjectWithCloudHandling :
  public Printable<GenericObjectWithCloudHandling> {
public:
  GenericObjectWithCloudHandling(bool Do_cloud = false)
    : do_cloud_(Do_cloud) {}
  virtual ~GenericObjectWithCloudHandling() { }

//-----------------------------------------------------------------------
/// If true, then do cloud retrieval, otherwise do clear.
//-----------------------------------------------------------------------

  bool do_cloud() const { return do_cloud_;}
  void do_cloud(bool F)
  {
    do_cloud_ = F;
    notify_do_cloud_update();
  }
  
//-----------------------------------------------------------------------
/// Derived class can override this to do an notifications needed whe
/// do_cloud_ is changed.
//-----------------------------------------------------------------------

  virtual void notify_do_cloud_update()
  {
  }
  virtual void print(std::ostream& Os) const 
  { Os << "GenericObjectWithCloudHandling"; }
protected:
  bool do_cloud_;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(GenericObjectWithCloudHandling);
#endif
