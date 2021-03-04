#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name="astro-sedpy",
    version='0.2.0',
    author="Ben Johnson",
    author_email="benjamin.johnson@cfa.harvard.edu",
    classifiers=["Development Status :: 4 - Beta",
                 "Intended Audience :: Science/Research",
                 "Programming Language :: Python",
                 "License :: OSI Approved :: MIT License",
                 "Natural Language :: English",
                 "Topic :: Scientific/Engineering :: Astronomy"],
    packages=["sedpy"],
    url="https://github.com/bd-j/sedpy",
    license="MIT",
    description=("Simple tools for astronomical spectral energy distributions, "
                 "particularly filter projections."),
    long_description=open("README.rst").read(),
    include_package_data=True,
    package_data={"sedpy": ["data/*fits", "data/filters/*par"]},
    install_requires=["numpy", "scipy", "astropy"],
)
