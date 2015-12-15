import os

import analysis
import bedtools
import biobambam
import cghub
import express
import featureCounts
import gencode
import general
import mutect
import picard
import pysamext
import star
import variants

_root = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]

scripts = os.path.join(
    os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-2]),
    'scripts')
