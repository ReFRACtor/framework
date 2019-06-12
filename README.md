This contains generic rules for wrapping code with SWIG.

This is not standalone code, it is meant to be included in a larger piece
of code. But this contains the SWIG specific part of this code.

Look for example at the [swig-rules-skeleton](https://github.jpl.nasa.gov/Cartography/swig-rules-skeleton) repository for a skeleton of using
this. A more full example can be found at 
[geocal](https://github.jpl.nasa.gov/Cartography/geocal).

Getting the SWIG rules into a repository
----------------------------------------

Normally you use git subtree to put this swig_rules. So assuming this has
been added as a remote swig-rules-repo you would do something like:

    git subtree pull --prefix swig_rules swig-rules-repo master
	
If you make local changes that should be shared with other repository, you
can go the opposite direction with:

    git subtree push --prefix swig_rules swig-rules-repo master
	
We can look at differences with upstream:

    git fetch swig-rules-repo
	git diff HEAD:swig_rules swig-rules-repo/master

The initial creation of the subtree in a new repository is by:

    git subtree add --prefix swig_rules swig-rules-repo master

Sometimes you will get errors saying "Updates were rejected because the
tip of current branch is behind." First do a subtree pull if you haven't,
this may be real. But you might get this even if a pull says everything
is up to date. This appears to be a bug in subtree, see description at
http://stackoverflow.com/questions/13756055/git-subtree-subtree-up-to-date-but-cant-push

Can try 

    git push swig-rules-repo `git subtree split --prefix swig_rules master`:master --force


Classes needed in code that uses this
-------------------------------------
The shared_ptr_type_mapper.i code depends on having a GenericObject 
available as generic_object.h. All the class should be derived from
this base class in order to be able to be cast automatically to the 
right type.

The code swig_type_mapper_base.cc should be compiled and included in
the parent library. This should have -DSWIG_MAPPER_NAMESPACE=GeoCal
or similar defined when compiled.
