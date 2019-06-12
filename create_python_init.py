# This is a short script for generating the top level wrapper code that
# initializes everything for SWIG/cpython. You don't normally run this directly,
# it is used by python.am for automatically creating this file.
#
# This takes a template, and a list of modules to combine. The list can contain
# things like "HAVE_HDF5", the modules listed after that will be included with
# conditional compilation "#ifdef HAVE_HDF5", so we can include modules based
# on various configuration options.
#
# First argument is the file to create, second is the template file to use.
#
# SWIG and cython are almost the same, so we pass "swig" or "cython" as the
# third argument.
#
# The fourth argument is the module name.
#
# The remaining arguments are the list of modules.

import os
import sys
import re

tmpl_dir = os.path.dirname(sys.argv[0]) + "/"
wrap_template = open(sys.argv[1]).read()
do_cython = (sys.argv[3] == "cython")
modname = sys.argv[4]
prototypes = []
initcmd = []
end_count = 0
for i in sys.argv[5:]:
    # Handle stuff we are including conditionally
    if(re.search('\AHAVE', i)):
        for c in range(end_count):
            prototypes.append("#endif")
            initcmd.append("#endif")
        end_count = 0
        for t in i.split():
            prototypes.append("#ifdef %s" % t)
            initcmd.append("#ifdef %s" % t)
            end_count += 1
    # Handle everything else
    else:
        if(do_cython):
            prototypes.append("  INIT_TYPE INIT_FUNC(%s)(void);" % i)
            initcmd .append("  INIT_MODULE(cython_list, \"%s.%s\", INIT_FUNC(%s));" % (modname, i, i))
        else:
            prototypes.append("  INIT_TYPE INIT_FUNC(_%s)(void);" % i)
            initcmd .append("  INIT_MODULE(module, \"_%s\", INIT_FUNC(_%s));" % (i, i))
# Make sure we close all the conditions
for c in range(end_count):
    prototypes.append("#endif")
    initcmd.append("#endif")
with open(sys.argv[2], 'w') as wrap_fo:
    wrap_fo.write(wrap_template.format(prototypes="\n".join(prototypes),
                                       initcmd="\n".join(initcmd)))


