# Check that all the dependent libraries are present. Should be the case,
# but we don't necessarily have all the libraries pinned. If this test fails,
# we should update the geocal recipe. In the mean time though, catch that
# an error occurred and give hints to the user about how to work around this

/usr/bin/env python $PREFIX/bin/.framework-post-link.py
