#!/usr/bin/env python
from obspy.core import read
import sys
import IPython
if __name__ == '__main__':
	file=sys.argv[1]
	aa=read(file,format='SAC')
	for kk in list(aa[0].stats.keys()):
		if not kk=='sac':
			try:
				print("\t"+kk+ ":"+ str(aa[0].stats.get(kk)))
			except Exception as ex1:
				IPython.embed()
				sys.exit()
	print('SAC header')
	for kk in aa[0].stats.get('sac'):
		try:
			if not float(aa[0].stats.get('sac').get(kk))==-12345:
				#IPython.embed()
				print("\t\t"+kk+':'+ str(aa[0].stats.get('sac').get(kk)))
		except Exception:
			print("\t\t"+kk+':'+ str(aa[0].stats.get('sac').get(kk)))
