import os
import bedtools, cghub, express, gencode, general, mutect, picard, pysamext
import star, variants

scripts = os.path.join(
    os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-2]),
    'scripts')
