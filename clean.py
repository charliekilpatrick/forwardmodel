#!/usr/bin/env python
import shutil
import sys
import os
import glob

def clean_dirs(basedir):
    for subdir in ['Medfilter','ReInterp','DualOutlier','BoxOutlier','Outlier',
        'Detect','Interp']:
        for dr in glob.glob(os.path.join(basedir, '*/sub_stack', subdir)):
            if os.path.exists(dr):
                shutil.rmtree(dr)

if __name__ == '__main__':
    if len(sys.argv)<2:
        basedir = sys.argv[1]
    else:
        print('Usage: clean.py basedir')
        sys.exit()

    if not os.path.exists(basedir):
        print(f'ERROR: basedir={basedir} does not exist')
        sys.exit()

    clean_dirs(basedir)


