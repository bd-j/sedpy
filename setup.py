#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

setup(
    name="sedpy",
    version='0.2.0',
    author="Ben Johnson",
    author_email="benjamin.johnson@cfa.harvard.edu",
    packages=["sedpy"],
    url="",
    license="LICENSE",
    description="Tools for dealing with astronomical spectral energy distributions",
    long_description=open("README.rst").read() + "\n\n",
    package_data={"sedpy": ["data/*fits", "data/filters/*par"]},
    include_package_data=True,
    install_requires=["numpy", "scipy >= 0.9", "astropy"],
)
