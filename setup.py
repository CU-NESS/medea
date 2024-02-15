#!/usr/bin/env python
"""
File: setup.py
Author: Joshua J. Hibbard
Date: Jan 2024

Description: Installs MEDEA.
"""
import os
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='medea', version='0.1',\
    description='Model to Emulate Directivities and Electric fields for Antennae',\
    author='Joshua J. Hibbard',\
    author_email='joshua.hibbard@colorado.edu',\
    packages=['medea'])

MEDEA_env = os.getenv('MEDEA')
cwd = os.getcwd()
if not MEDEA_env:
    import re
    shell = os.getenv('SHELL')
    print("\n")
    print("#" * 78)
    print("It would be in your best interest to set an environment variable")
    print("pointing to this directory.\n")
    if shell:
        if re.search('bash', shell):
            print("Looks like you're using bash, so add the following to " +\
                "your .bashrc:")
            print("\n    export MEDEA={0}".format(cwd))
        elif re.search('csh', shell):
            print("Looks like you're using csh, so add the following to " +\
                "your .cshrc:")
            print("\n    setenv MEDEA {!s}".format(cwd))
    print("\nGood luck!")
    print("#" * 78)
    print("\n")
elif MEDEA_env != cwd:
    print("\n")
    print("#" * 78)
    print("It looks like you've already got an medea environment variable " +\
        "set but it's \npointing to a different directory:")
    print("\n    MEDEA={!s}".format(MEDEA_env))
    print("\nHowever, we're currently in {!s}.\n".format(cwd))
    print("Is this a different medea install (might not cause problems), or " +\
        "perhaps just")
    print("a typo in your environment variable?")
    print("#" * 78)
    print("\n")
