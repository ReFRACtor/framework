The lifetime of python classes that derive from C++ director objects is a bit
confusing, since we have both python and C++ pieces it isn't immediately
clear who owns things, and how lifetime/cleanup is handled.

After spending some time with this, we have adopted python *always* owns the
objects and manages the lifetimes. A python director object is given a 
"this" member function by SWIG when we create the director object. The
"this" member points to a boost::shared_ptr of the underlying C++ code
(so in our skeleton [swig-rules-skeleton](https://github.jpl.nasa.gov/Cartography/swig-rules-skeleton) DirectorExample with is a C++ class 
SwigDirector_DirectorExample object).

When the python object gets a zero reference count, it is deleted and it in
turns deletes the C++ SWIG director object.

The complication that arises is when we pass this python object to C++ where
it appears as the same SwigDirector_DirectorExample class, which contains
a pointer to the python class (as swig_self member variable, available by
the C++ accessor code swig_get_self()), which in turns contains the python
level "this" that points back to the same SwigDirector_DirectorExample.

The way we handle this is that we handle the passed C++ class in SWIG as
a *separate* shared_ptr - so a different shared_ptr than the one used in
the "this" python object. This shared_ptr used a feature of std::shared_ptr
where you can supply a deleter object. The default shared_ptr deletes the
data it points to when its count gets to zero, but this can behavior can
be changed by supplying the deleter. We use PythonRefPtrCleanup class, which
instead of deleting the SwigDirector_DirectorExample object instead decrements
the python object reference count.

If we happen to only have one reference left to the python director class
as a shared_ptr to the C++ SWIG object, then when the SWIG object isn't
directly deleted. Instead, the python reference is decremented, results in
a count of 0, which then results in the python object cleaning of "this"
and finally deleting the actual C++ SWIG object.

We modify both the swig code (so my_shared_ptr.i) and the boost serialization
code of directors to handle this.

See the [swig-rules-skeleton](https://github.jpl.nasa.gov/Cartography/swig-rules-skeleton) for an example of this in use, and in particular py-tests to
check that the object lifetimes are handled correctly.

Handling weak_ptr
-----------------

There is a corner case that needs to be handled. Since we normally return a *different* 
shared_ptr, this breaks the handling of weak_ptr. Without special handling, a
shared_ptr is passed to a function setting a weak_ptr. This shared_ptr has a ref_count
of 1. The weak_ptr is set to point to this. The function ends, and the shared_ptr has
its ref_count set to 0. The underlying object still existing, PythonRefPtrCleanup handles
the lifetime and python still have references to it so the object isn't deleted. However,
the C++ weak_ptr sees a shared_ptr with a 0 reference count, so it drops the reference to
the object. We have a unit test illustrating this problem in swig-rules-skeleton, see
DirectorExampleWeakPtr.

So we have a special typemap that triggers off the variable name SHARED_PTR_NO_OWN. This 
returns the shared_ptr owned by the python object for Director. This shared pointer should
not be used for anything we keep a copy of because this will break the lifetime issues. But
it can be used to set a weak_ptr. For a normal C++ object that isn't a python object
Director this is handled the same as the normal typemap rule for shared_ptr.

Note that it is important that SHARED_PTR_NO_OWN not actually be copied to a persistent 
shared_ptr. If you do, then you are likely to get a seg fault - this is one of those
"giving you enough rope to hang yourself" sort of thing. This handles a corner case, but
by *breaking* in a controlled way the lifetime management that we have in place.

BTW, we orignally put the weak_ptr handling in place for being able to
have ReFRACtor framework Observer working correctly. Turns out this
wasn't actually needed there. We *already* handled this a different
way in Observer, I had forgotten this. In that case the handling was
accidental, we wanted to be able to have Observer passed in as
something other than a shared_ptr, so we added a "this_obj" in it that
had the same lifetime as the containing object - and made a weak_ptr
to this. So although we have this fixed in our swig rules now with a
special typemap, it turns out we didn't actually need this for
ReFRACtor. Might be useful in some other context though, so it was
still a worthwhile excercise.

