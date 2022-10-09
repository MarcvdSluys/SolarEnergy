#!/bin/env python
# -*- coding: utf-8 -*-

"""Setup.py for the SolarEnergy Python package."""

# Package version:
version='0.1.3'

# Get long description from README.md:
with open('README.md', 'r') as fh:
    long_description = fh.read()


from setuptools import setup
setup(
    name='solarenergy',
    description='A Python module do simple modelling in the field of solar energy',
    author='Marc van der Sluys',
    url='https://github.com/MarcvdSluys/SolarEnergy',
    
    packages=['solarenergy'],
    install_requires=['astroconst>=0.0.5','numpy','pytz','soltrack>=0.2.0'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    version=version,
    license='EUPL 1.2',
    keywords=['solar energy','solar','energy','sun'],
    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
    ]
)
