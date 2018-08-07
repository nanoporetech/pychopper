#!/usr/bin/env python

"""Helper script to generate cmd_line.rst file for all scripts in bin which
have and parser object defined in their global scope - taken from tang

"""
from __future__ import print_function
import sys
import os
import imp

scripts_rel = 'scripts'
attr_name = 'parser'
blacklist = ['__init__.py']

location = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
scripts_abs = os.path.join(location, scripts_rel)
scripts = sorted(filter(lambda s: s[0] != '.' and s not in blacklist, os.listdir(scripts_abs)))

sys.stderr.write("Found following scripts:\n{}\n{}\n{}\n".format(
    location, scripts_abs, scripts
))


print ("""
.. _command_line_tools:

Command line tools
==================
""")

for script in scripts:
    script_name, script_ext = os.path.splitext(script)
    if script_ext == '.pyc':
        continue
    
    try:
        mod_name = '{}.{}'.format(scripts_rel, script_name)
        #mod = __import__(mod_name, globals(), locals(), [attr_name])
        mod = imp.load_source(script_name, os.path.join(scripts_abs, script))
        script = script.replace('.py', '')

        print ('.. _{}:\n\n{}\n{}'.format(script, script, '-'*len(script)))
        if hasattr(mod, attr_name):
            print ("""
.. argparse::
   :ref: {}.{}
   :prog: {}
""".format(mod_name, attr_name, script_name))
        else:
            print('No documentation available')

    except Exception as e:
        # Wha' yer' gonna do?
        sys.stderr.write('Error making docs for {}:\n{}\n'.format(script_name, e))
        pass
