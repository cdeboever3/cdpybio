import os

__version__ = '0.0.8'

__all__ = [
    "analysis",
    "bedtools",
    "biobambam",
    "express",
    "featureCounts",
    "gencode",
    "general",
    "ldsc",
    "moodsext",
    "mutect",
    "picard",
    "plink",
    "pysamext",
    "star",
    "variants",
]

from . import analysis
from . import bedtools
from . import biobambam
from . import express
from . import featureCounts
from . import gencode
from . import general
from . import ldsc
from . import moodsext
from . import mutect
from . import picard
from . import plink
from . import pysamext
from . import star
from . import variants

_root = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]

scripts = os.path.join(
    os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-2]),
    'scripts')
