import os
import sys
from setuptools import setup, find_packages
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/src')
from src import __version__
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
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    packages=find_packages(),
    install_requires = [
        'pysam>=0.14.0',
        'pybedtools>=0.9.0',
        'numpy>=1.22.3',
        'pathos>=0.2.3',
        'scikit-learn>=1.1',
        'scipy>=1.8.0',
        ],
    scripts = ['straglr.py',
               ],
)
