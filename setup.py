import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst', format='md')
except(IOError, ImportError):
    long_description = read('README.md')

setup(
    name = 'cdpybio',
    packages=['cdpybio'],
    version = '0.0.4',
    author = 'Christopher DeBoever',
    author_email = 'cdeboever3@gmail.com',
    description = ('A collection of useful methods for working with various '
                   'bioinformatics data, software output files, etc.'),
    license = 'MIT',
    keywords = ['bioinformatics'],
    url = 'https://github.com/cdeboever3/cdpybio',
    download_url = 'https://github.com/cdeboever3/cdpybio/tarball/0.0.4',
    long_description=long_description,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
   ]
)
