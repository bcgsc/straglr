import os
from setuptools import setup, find_packages
from straglr import __version__
setup(
    name='straglr',
    version=__version__,
    description='Straglr',
    long_description='Short tandem repeat genotyping using long reads',
    url='https://github.com/bcgsc/straglr.git',
    author='Readman Chiu',
    author_email='rchiu@bcgsc.ca',
    license='BCCA',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    packages=find_packages(),
    install_requires = [
        'pysam>=0.14.0',
        'pybedtools>=0.7.9',
        'intspan>=1.5.8',
        'numpy>=1.14.2',
        'pathos>=0.2.3',
        'scikit-learn>=0.19.0',
        'scipy>=1.0.1',
        ],
    scripts = ['straglr.py',
               ],
)
