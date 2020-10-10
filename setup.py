#!/bin/env python3

"""Setup.py for the SolarEnergy Python package."""


# Package version:
version="0.0.12"

# Get long description from README.md:
with open("README.md", "r") as fh:
    long_description = fh.read()


from setuptools import setup
setup(
    name='solarenergy',
    description='A Python module do simple modelling in the field of solar energy',
    author='Marc van der Sluys',
    url='https://github.com/MarcvdSluys/SolarEnergy',
    
    packages=['solarenergy'],
    install_requires=['pytz','numpy','soltrack'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    version=version,
    license='GPLv3+',
    keywords=['solar energy','solar','energy','sun'],
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
    ]
)

