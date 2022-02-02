from setuptools import setup, find_packages
import pathlib

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst', format='md')
except(IOError, ImportError):
    long_description = open('README.md').read()


here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

install_requires = [
    "biopython>=1.72",  # BioPython
    "HTSeq>=0.8.0",
    "numpy>=1.15.0",
    "scipy>=1.4.1",
    "pandas>0.23.4",
    "pysam>=0.15.1",
    "pybedtools>0.8.0",
    "PyVCF>=0.6.8",
    "pytest>=5.0.0",
    "seaborn>0.9.0",
    "rpy2>=3.0.0",
    # "pdb", #  ??
    "statsmodels>=0.8.0",
    "matplotlib>3.0.0",
    "gffutils>=0.9",
    "urllib3>=1.25.0",
    # "MOODS>=1.9.4.1", #  https://github.com/jhkorhonen/MOODS (2019) Python
]

setup(
    name = 'cdpybio',
    packages=['cdpybio'],
    version = '0.0.8',
    author = 'Christopher DeBoever',
    author_email = 'cdeboever3@gmail.com',
    description = ('A collection of useful methods for working with various '
                   'bioinformatics data, software output files, etc.'),
    license = 'MIT',
    keywords = ['bioinformatics'],
    url = 'https://github.com/FredHutch/cdpybio',
    #download_url = 'https://github.com/cdeboever3/cdpybio/tarball/0.0.7',
    long_description_content_type='text/markdown',
    long_description=long_description,
    install_requires=install_requires,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        "Programming Language :: Python :: 3.10",
        'Programming Language :: Python :: 3 :: Only',
    ],
    python_rquires='>=3.6',
    project_urls={
        'Bug Reports': 'https://github.com/FredHutch/cdpybio/issues',
        'Source': 'https://github.com/FredHutch/cdpybio',
    },
)
