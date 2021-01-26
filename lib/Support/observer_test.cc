#include "observer.h"
#include "unit_test_support.h"

using namespace FullPhysics;

class ObservableTest : public Observable<ObservableTest> {
public:
  ObservableTest() :v_(0) {}
  virtual ~ObservableTest() {}
  int v() const { return v_;}
  void v(int newv) { v_ = newv; notify_update_do(*this);}
  virtual void add_observer(Observer<ObservableTest>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<ObservableTest>& Obs) 
  { remove_observer_do(Obs, *this);}
private:
  int v_;
};

class CacheInvalidatedObserverTest : public CacheInvalidatedObserver {
public:
  CacheInvalidatedObserverTest() {}
  virtual ~CacheInvalidatedObserverTest() {}
  bool& cachev() { return cache_valid; }
};
  
class ObserverTest : public Observer<ObservableTest> {
public:
  ObserverTest() : myv_(-1) {}
  virtual ~ObserverTest() {}
  int myv() { return myv_; }
  virtual void notify_update(const ObservableTest& T)
  {
    myv_ = T.v();
  }
private:
  int myv_;
};

BOOST_FIXTURE_TEST_SUITE(observer, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  ObservableTest t;
  ObserverTest t2;
  CacheInvalidatedObserverTest t3;
  t3.cachev() = true;
  t.add_observer(t2);
  t.add_cache_invalidated_observer(t3);
  BOOST_CHECK_EQUAL(t2.myv(), -1);
  BOOST_CHECK_EQUAL(t3.cachev(), true);
  t.v(10);
  BOOST_CHECK_EQUAL(t2.myv(), 10);
  BOOST_CHECK_EQUAL(t3.cachev(), false);
}

BOOST_AUTO_TEST_CASE(proper_cleanup)
{
  ObservableTest t;
  boost::shared_ptr<ObserverTest> t2(new ObserverTest);
  boost::shared_ptr<CacheInvalidatedObserverTest>
    t3(new CacheInvalidatedObserverTest);
  t3->cachev() = true;
  t.add_observer(*t2);
  t.add_cache_invalidated_observer(*t3);
  BOOST_CHECK_EQUAL(t2->myv(), -1);
  BOOST_CHECK_EQUAL(t3->cachev(), true);
  t.v(10);
  BOOST_CHECK_EQUAL(t2->myv(), 10);
  BOOST_CHECK_EQUAL(t3->cachev(), false);
  t2.reset();
  t3.reset();
  t.v(10);
}

BOOST_AUTO_TEST_SUITE_END()
