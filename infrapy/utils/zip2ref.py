#!/usr/bin/env python

import glob
import sys
import subprocess
from  IPython import embed
directory=sys.argv[1]

files=glob.glob(directory+'/*.zip')
print(files)

for ii in files:
    print(ii)
    iin=ii.split('/')
    donefile=glob.glob('zip_conv/'+iin[-1]+'done')
    #embed()
    #exit()
    if len(donefile)==1:
        print(ii+ '  already converted')
        continue
    #embed()
    #exit()
    try:
        subprocess.check_output('tar -xkf \''+ii+'\' --strip-components=1 -C ref_files/;',shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        #embed()
        #exit()
    try:
        subprocess.check_output('touch zip_conv/\''+iin[-1]+'\'done;',shell=True)
        except subprocess.CalledProcessError as e:
                exit()
