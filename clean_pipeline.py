#!/usr/bin/env python
import shutil
import sys
import os
import glob

def clean_dirs(basedir):
    for file_glob in ['*.pickle', '*.txt', '*.log', '*.sh', '*.fits']:
        for file in glob.glob(os.path.join(basedir, file_glob)):
            if not os.path.isfile(file): continue
            print(f'Deleting file: {file}')
            os.remove(file)

    for subdir_glob in ['sub_stack_*','psfs','*:*','r*','ut*']:
        for subdir in glob.glob(os.path.join(basedir, subdir_glob)):
            if not os.path.isdir(subdir): continue
            print(f'Deleting directory: {subdir}')
            shutil.rmtree(subdir)


if __name__ == '__main__':
    if len(sys.argv)>1:
        basedir = sys.argv[1]
    else:
        basedir = '.'

    if not os.path.exists(basedir):
        print(f'ERROR: basedir={basedir} does not exist')
        sys.exit()

    clean_dirs(basedir)


