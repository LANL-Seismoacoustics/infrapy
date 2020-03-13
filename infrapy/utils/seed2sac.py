from obspy.clients.iris import Client
from obspy import UTCDateTime
from obspy.io.xseed import Parser
from obspy.core import read

from IPython import embed

from obspy.core.util import AttribDict
import argparse


# get seed from

import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':
    par = argparse.ArgumentParser(description="Convert a full seed volume to sac")
    par.add_argument('-d', dest='seed',required=True,help='seed volume, e.g.: -d tes1.seed')
    par.add_argument('-l', dest='length',required=False,help='length in seconds of individual sac files, e.g.: -l 7200')
    args=par.parse_args()
    if args.seed:
        seed=args.seed
    if args.length:
        length=float(args.length)
    else:
        print('using 7200 s for length of sac files')
        length=7200
    p = Parser(seed)
    aa =read(seed)
    for ai in aa:
        print(ai)
        st=ai.stats.starttime
        et=ai.stats.endtime
        ct=st
        while ct < et:
            ai_aux=ai.copy()
            ai_aux.trim(starttime=ct,endtime=ct+length)
            pos=p.get_coordinates(ai_aux.get_id())
            ai_aux.stats.sac=AttribDict({'stla':pos['latitude'],'stlo': pos['longitude'],'stel': 0})
            ai_aux.write(ai_aux.stats.station+str(ct)+'.sac',format='SAC')
            ct=ct+length
            print(str(ct))
