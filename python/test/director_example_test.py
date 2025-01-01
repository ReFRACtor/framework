from refractor.framework import *
from test_support import *
import pickle
import sys
import gc
import ctypes

# -----------------------------------------------------------------------------------
# Note, this is a duplicate of a test we have in the much smaller swig-rules-skeleton.
# Unfortunately the handling of director objects by swig doesn't fully handle shared
# pointers. We have fixes for this in our rules, but this occasionally breaks as the
# actual swig implementation changes. We've had this break 3 times now as we've moved
# to newer swig versions.
#
# It seems desirable to fail quickly - the symptom of the director failures is that
# we suddenly have seg faults in various places. So it would be good to have a unit
# test that fails rather than catching this later when using framework.
#
# If you get a failure here, you can debug this in framework but it is generally easier
# to work in swig-rules-skeleton which is a much smaller system that just tests the
# swig interface.
#
# You can usually work around a problem by just temporarily requiring an older version
# of swig until the underlying problem can be solved. But you generally need to debug
# the low level swig wrappers written by swig, so fixing this isn't trivial. But we
# have a lot of support in swig-rules-skeleton for fixing this if needed.
# -----------------------------------------------------------------------------------


# For testing, this calls the low level python c API for decreasing the reference
# count. Can use to try to sort out issues with lifetime
_decref = ctypes.pythonapi.Py_DecRef
_decref.argtypes = [ctypes.py_object]
_decref.restype = None


class DirectorPythonExample(rf.DirectorExample):
    # Note, tests that use deleted_list need to run in serial
    deleted_list = set()

    def __init__(self, v1, v2, name):
        super().__init__(v1, name)
        self.v2 = v2
        print(f"Creating Python {self.name()}", flush=True)

    def __del__(self):
        self.deleted_list.add(id(self))
        print(f"Destroying Python {self.name()}", flush=True)

    def print_ref_count(self):
        # We have an extra 2 counts because this got passed as an
        # arg to getrefcount, and we got this as a arg.
        print(f"Ref count {self.name()} {sys.getrefcount(self)-2}")
        
    def func(self, v):
        return super().func(v) + 2 * self.v2

    def desc(self):
        '''Description to print for this object.'''
        return f'''DirectorPythonExample:
   name:                { self.name() }
   ref_count:           { sys.getrefcount(self)-2 }
   value_to_add:        { self.value_to_add() }
   value_to_double_add: {self.v2}'''

    def print_obj():
        '''Print a list of active DirectorPythonExample, and the first one if it
        exists. Useful when sorting through lifetimes issues.'''
        olist = [obj for obj in gc.get_objects() if isinstance(obj, DirectorPythonExample)]
        print(f"Found in active objects: {olist}", flush=True)
        if(len(olist) > 0):
            print("First object in active objects:")
            if(olist[0].this is not None):
                print(olist[0])
                print(gc.get_referents(olist[0].this), flush=True)
                print(sys.getrefcount(olist[0].this))
            else:
                print(sys.getrefcount(olist[0]))
            print("Stuff referring to this:")
            print(gc.get_referents(olist[0]), flush=True)
    
def test_director_example(serial):
    print("May want to valgrind this, to make sure we have lifetime etc. ok with memory usage.")
    print("Do this with PYTHONMALLOC=malloc valgrind $(which python) $(which pytest) -s lib/director_example_test.py")
    print("(Note PYTHONMALLOC is important, otherwise you get a ton of spurious errors)")
    DirectorPythonExample.deleted_list.clear()
    e1 = rf.DirectorExample(5, "e1")
    assert e1.func(10) == 5+10
    e1 = None
    e2 = DirectorPythonExample(5,6, "e2")
    assert e2.func(10) == 5+10+2*6
    euser = rf.DirectorExampleUser()
    euser.set_value(e2)
    e2 = None
    assert euser.apply_func(10) == 5+10+2*6
    print(euser.value())
    # Test a simple std::vector
    euser.set_vec_int([1,2,3])
    print(euser.value_vec_int())
    if True:
        # Test a std::vector of boost::shared_ptr
        v = [rf.DirectorExample(5, "v1"),
             rf.DirectorExample(6, "v2"),
             DirectorPythonExample(5,6, "vp3"),
             DirectorPythonExample(10,7, "vp4")]
    else:
        # Directly use std::vector. Shouldn't be needed any longer, passing
        # in a list works now. But leave this around in case we run into a
        # future issue we need to diagnose
        v = rf.Vector_DirectorExample()
        v.push_back(rf.DirectorExample(5, "v1"))
        v.push_back(rf.DirectorExample(6, "v2"))
        v.push_back(DirectorPythonExample(5,6, "vp3"))
        v.push_back(DirectorPythonExample(10,7, "vp4"))
    euser.set_vec(v)
    v = None
    assert euser.apply_func(10, 0) == 5+10
    assert euser.apply_func(10, 1) == 6+10
    assert euser.apply_func(10, 2) == 5+10+2*6
    assert euser.apply_func(10, 3) == 10+10+2*7
    print(euser.value_vec()[0])
    print(euser.value_vec()[1])
    print(euser.value_vec()[2])
    print(euser.value_vec()[3])
    assert euser.value_vec()[0].__class__ == rf.DirectorExample
    assert euser.value_vec()[1].__class__ == rf.DirectorExample
    assert euser.value_vec()[2].__class__ == DirectorPythonExample
    assert euser.value_vec()[3].__class__ == DirectorPythonExample
    print("Should see e2 and v1 through vp4 deleted", flush=True)
    assert len(DirectorPythonExample.deleted_list) == 0
    euser = None
    assert len(DirectorPythonExample.deleted_list) == 3
    DirectorPythonExample.deleted_list.clear()
    print("All Done", flush=True)
    euser = rf.DirectorExampleUser()
    v2 = [[rf.DirectorExample(5, "v1"), DirectorPythonExample(5,6, "vp2")],
          [rf.DirectorExample(6, "v2"), DirectorPythonExample(10,7, "vp4"),
           DirectorPythonExample(11,7, "vp5")]]
    euser.set_vec_vec(v2)
    v2 = None
    print(euser.value_vec_vec())
    assert euser.value_vec_vec()[0][0].__class__ == rf.DirectorExample
    assert euser.value_vec_vec()[0][1].__class__ == DirectorPythonExample
    assert euser.value_vec_vec()[1][0].__class__ == rf.DirectorExample
    assert euser.value_vec_vec()[1][1].__class__ == DirectorPythonExample
    assert euser.value_vec_vec()[1][2].__class__ == DirectorPythonExample
    assert euser.apply_func(10, 0, 0) == 5+10
    assert euser.apply_func(10, 1, 0) == 6+10
    assert euser.apply_func(10, 0, 1) == 5+10+2*6
    assert euser.apply_func(10, 1, 1) == 10+10+2*7
    assert euser.apply_func(10, 1, 2) == 10+11+2*7
    print("Should see v1 through vp5 deleted", flush=True)
    assert len(DirectorPythonExample.deleted_list) == 0
    euser = None
    assert len(DirectorPythonExample.deleted_list) == 3
    DirectorPythonExample.deleted_list.clear()
    print("All Done", flush=True)

def test_director_pickle_shelve(isolated_dir, serial):
    print("May want to valgrind this, to make sure we have lifetime etc. ok with memory usage.")
    print("Do this with PYTHONMALLOC=malloc valgrind $(which python) $(which pytest) -s lib/director_example_test.py")
    print("(Note PYTHONMALLOC is important, otherwise you get a ton of spurious errors)")
    e1 = rf.DirectorExample(5, "e1")
    assert e1.func(10) == 5+10
    # This is tiny, but go ahead and compress just to test the interface
    write_shelve("e1.xml.gz", e1)
    e1_s = read_shelve("e1.xml.gz")
    pickle.dump(e1, open("e1.pkl", "wb"))
    e1_p = pickle.load(open("e1.pkl", "rb"))
    e1 = None
    assert e1_s.func(10) == 5+10
    assert e1_p.func(10) == 5+10
    print(e1_s)
    print(e1_p)
    e1_s = None
    e1_p = None
    e2 = DirectorPythonExample(5,6, "e2")
    assert e2.func(10) == 5+10+2*6
    write_shelve("e2.xml", e2)
    e2_s = read_shelve("e2.xml")
    pickle.dump(e2, open("e2.pkl", "wb"))
    e2_p = pickle.load(open("e2.pkl", "rb"))
    assert len(DirectorPythonExample.deleted_list) == 0
    e2 = None
    assert len(DirectorPythonExample.deleted_list) == 1
    DirectorPythonExample.deleted_list.clear()
    assert e2_s.func(10) == 5+10+2*6
    assert e2_p.func(10) == 5+10+2*6
    print(e2_s)
    print(e2_p)
    assert len(DirectorPythonExample.deleted_list) == 0
    e2_s = None
    e2_p = None
    assert len(DirectorPythonExample.deleted_list) == 2
    DirectorPythonExample.deleted_list.clear()
    
    euser = rf.DirectorExampleUser()
    v = [rf.DirectorExample(5, "v1"),
         rf.DirectorExample(6, "v2"),
         DirectorPythonExample(5,6, "vp3"),
         DirectorPythonExample(10,7, "vp4")]
    euser.set_vec(v)
    v = None
    write_shelve("euser.xml", euser)
    euser_s = read_shelve("euser.xml")
    pickle.dump(euser, open("euser.pkl", "wb"))
    euser_p = pickle.load(open("euser.pkl", "rb"))
    assert len(DirectorPythonExample.deleted_list) == 0
    euser = None
    assert len(DirectorPythonExample.deleted_list) == 2
    DirectorPythonExample.deleted_list.clear()
    assert euser_s.apply_func(10, 0) == 5+10
    assert euser_s.apply_func(10, 1) == 6+10
    assert euser_s.apply_func(10, 2) == 5+10+2*6
    assert euser_s.apply_func(10, 3) == 10+10+2*7
    assert euser_p.apply_func(10, 0) == 5+10
    assert euser_p.apply_func(10, 1) == 6+10
    assert euser_p.apply_func(10, 2) == 5+10+2*6
    assert euser_p.apply_func(10, 3) == 10+10+2*7
    print(euser_s.value_vec()[0])
    print(euser_s.value_vec()[1])
    print(euser_s.value_vec()[2])
    print(euser_s.value_vec()[3])
    print(euser_p.value_vec()[0])
    print(euser_p.value_vec()[1])
    print(euser_p.value_vec()[2])
    print(euser_p.value_vec()[3])
    assert len(DirectorPythonExample.deleted_list) == 0
    euser_s = None
    euser_p = None
    assert len(DirectorPythonExample.deleted_list) == 4
    DirectorPythonExample.deleted_list.clear()
    print("All done", flush=True)

def test_director_lifetime_shelve(isolated_dir, serial):
    DirectorPythonExample.deleted_list.clear()
    print("Should see creation message", flush=True)
    e1 = DirectorPythonExample(5,6, "e1")
    e1_id = id(e1)
    write_shelve("e1.xml", e1)
    print("Should see deletion message", flush=True)
    assert e1_id not in DirectorPythonExample.deleted_list
    e1 = None
    assert e1_id in DirectorPythonExample.deleted_list
    DirectorPythonExample.deleted_list.clear()
    print("All done with e1", flush=True)
    e1_s = read_shelve("e1.xml")
    e1_s_id = id(e1_s)
    print(e1_s)
    print("Should see deletion message", flush=True)
    assert e1_s_id not in DirectorPythonExample.deleted_list
    e1_s = None
    assert e1_s_id in DirectorPythonExample.deleted_list
    DirectorPythonExample.deleted_list.clear()
    print("All done with e1_s", flush=True)
    print("Should see creation message", flush=True)
    euser = rf.DirectorExampleUser()
    euser.set_vec([DirectorPythonExample(5,6, "e2"),])
    # We can't easily get the id DirectorPythonExample, because it
    # gets copied when we return it. So just count the number of
    # items in the deletion set, that should be pretty accurate 
    #e2_vec_id = id(euser.value_vec()[0])
    write_shelve("euser.xml", euser)
    print("Should see deletion message", flush=True)
    assert len(DirectorPythonExample.deleted_list) == 0
    euser = None
    assert len(DirectorPythonExample.deleted_list) == 1
    DirectorPythonExample.deleted_list.clear()
    print("All done with euser")
    euser_s = read_shelve("euser.xml")
    print(euser_s.value_vec()[0])
    print("Should see deletion message", flush=True)
    assert len(DirectorPythonExample.deleted_list) == 0
    euser_s = None
    assert len(DirectorPythonExample.deleted_list) == 1
    DirectorPythonExample.deleted_list.clear()
    print("All done with euser_s")

@skip    
def test_director_example_weak_ptr(isolated_dir, serial):
    m = GenericObjectMap()
    m["e1"] = DirectorPythonExample(5,6, "e1")
    m["e1_weak"] = DirectorExampleWeakPtr(m.e1)
    print(m.e1)
    print(m.e1_weak)
    assert m.e1.func(10) == 5+10+2*6
    assert m.e1_weak.a().func(10) == 5+10+2*6
    write_shelve("e1.xml", m)
    m_s = read_shelve("e1.xml")
    assert m_s.e1.func(10) == 5+10+2*6
    assert m_s.e1_weak.a().func(10) == 5+10+2*6
    
    # Check that m_s.b.a and m_s.a are the same. We can't directly compare, python
    # has two different shared ptr so the are "different". But we can verify
    # that the underlying C++ object is the same by changing one and making
    # sure the other is changed also
    m_s.e1.value_to_add(3)
    assert m_s.e1.func(10) == 3+10+2*6
    assert m_s.e1_weak.a().func(10) == 3+10+2*6
    m_s.e1.v2 = 7
    assert m_s.e1.func(10) == 3+10+2*7
    assert m_s.e1_weak.a().func(10) == 3+10+2*7

    # Do the same, but without using a GenericObjectMap. In this case, 
    # the weak_ptr in C should be set to None.
    a = DirectorPythonExample(5,6, "e1")
    c = rf.DirectorExampleWeakPtr(a)
    write_shelve("a.xml", a)
    a2 = read_shelve("a.xml")
    write_shelve("c.xml", c)
    c2 = read_shelve("c.xml")
    assert a2.func(10) == 5+10+2*6
    assert c2.a() is None

    # Also check handling of weak_ptr
    assert c.a() is not None
    assert not c.expired()
    assert c.use_count() == 1
    a = None
    assert c.expired()
    assert c.use_count() == 0
    assert c.a() is None
    

    
