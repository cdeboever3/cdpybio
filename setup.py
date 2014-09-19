import os
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def post(command_subclass):
    """A decorator for classes subclassing one of the setuptools commands.

    It modifies the run() method to download firehose_get
    """
    orig_run = command_subclass.run

    def modified_run(self):
        import shutil
        from subprocess import check_call
        from urllib2 import urlopen
        src = 'http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip'
        req = urlopen(src)
        dest = 'bin/firehose_get_latest.zip'
        try:
            os.makedirs('bin')
        except OSError:
            pass
        with open(dest, 'w') as d:
            shutil.copyfileobj(req, d)
        check_call(['unzip', dest, '-d', 'bin'])
        os.remove(dest)
        orig_run(self)

    command_subclass.run = modified_run
    return command_subclass

@post
class CustomDevelopCommand(develop):
    pass

@post
class CustomInstallCommand(install):
    pass

setup(
    name = 'cdpybio',
    version = '0.1.0',
    author = 'Christopher DeBoever',
    author_email = 'cdeboever3@gmail.com',
    description = ('A module that contains different useful things I like to '
                   'use.'),
    packages=find_packages(),
    cmdclass={
        'install' : CustomInstallCommand,
        'develop' : CustomDevelopCommand,
    },
    license = 'MIT',
    keywords = 'bioinformatics',
    url = 'https://github.com/cdeboever3/cdpybio',
    long_description=read('README.md'),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
   ]
)
