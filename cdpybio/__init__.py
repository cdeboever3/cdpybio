import os
import analysis, bedtools, biobambam, cghub, express, gencode, general, mutect
import picard, pysamext, star, variants

_root = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]

scripts = os.path.join(
    os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-2]),
    'scripts')
