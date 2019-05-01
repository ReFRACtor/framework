import subprocess
import os
import re
import sys

# For now, only run on linux
if(sys.platform == "darwin"):
    exit(0)
prefix = os.environ["PREFIX"]
libname = "%s/lib/librefractor.so" % prefix
t = subprocess.run("ldd %s | grep 'not found' | sort | uniq " % libname,
                   shell=True, stdout=subprocess.PIPE)
if(t.stdout == b''):
    exit(0)
fh = open("%s/.messages.txt" % prefix, 'a')
print("The following libraries are missing:\n", file=fh)  
for missing in t.stdout.split(b'\n'):
    m = re.match(rb'\s*(.*) => not found\s*', missing)
    if(m):
        print("  %s" % m.group(1).decode('utf-8'), file=fh)
print('''
Conda doesn't do a great job of tracking library dependencies. It
seems to assume that library versions only change at major package
version changes, i.e., version 11.0 and 11.2 of a package have the
same library version, but 12.0 is different. Many packages violate this
assumption. Conda can be explicitly told in a package build recipe to pin the
dependent package versions to a certain range, but this needs to be
done explicitly.

We've tried to pin packages we know are a problem in ReFRACtor
Framework, but it is an ongoing effort. You don't want to over specify
the package versions, because that can lead to breakage with other
packages (e.g., framework depends on foo version 1.11, but
awesome_package depends on foo version 1.12).

We check on installing ReFRACtor Framework if all the dependency
libraries are found.  We've failed this check, which is why you are
seeing this message.

This can often be worked around by explicitly installing a version
of the dependent package that works with framework in an environment
before installing framework (e.g., "conda create -n framework-test framework
hdf5=1.10.4").

Make sure to send a bug report for framework. You shouldn't need to
specify the version of dependency to install. The whole point of a
packaging system like conda is to track those dependencies for
you. But giving an explicit version will let you work around the
problem until the conda build recipe is improved.  
''', file=fh)
exit(1)


