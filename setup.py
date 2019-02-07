#!/usr/bin/env python
# coding=utf-8
# This file is part of pySOS.
#
# Copyright 2019 Francisco Nemiña

import os
from setuptools import setup


PROJECT_ROOT = os.path.dirname(__file__)


def read_file(filepath, root=PROJECT_ROOT):
    """
    Return the contents of the specified `filepath`.
    * `root` is the base path and it defaults to the `PROJECT_ROOT` directory.
    * `filepath` should be a relative path, starting from `root`.
    """
    with open(os.path.join(root, filepath)) as fd:
        text = fd.read()
    return text


LONG_DESCRIPTION = read_file("README.md")
SHORT_DESCRIPTION = "pySOS is a python interface for the Successive Orders of"\
                    + "Scattering (SOS)"\
                    + "radiative transfer code."
REQS = [
    'numpy',
    'matplotlib',
    'scipy'
]


setup(
    name                  = "pySOS",
    packages              = ['pySOS'],
    install_requires      = REQS,
    version               = "0.1a",
    author                = "Francisco Nemiña",
    author_email          = "fnemina@conae.gov.ar",
    description           = SHORT_DESCRIPTION,
    license               = "GPL-3.0",
    url                   = "None",
    long_description      = LONG_DESCRIPTION,
    classifiers           = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 2"
    ],
)
