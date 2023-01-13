import yaml
import IODDC10 as iod
import numpy
import subprocess, os
from datetime import datetime
import argparse

parser = argparse.ArgumentParser(description='wrapper for DDC10 running',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('conf',type=str,help="Path to yaml config file for DDC10")

args=parser.parse_args()

mDDC10 = iod.IODDC10.from_yml(args.conf)
rPar = mDDC10.config
mDDC10.setupDDC10(fade=rPar['fade'],force=rPar['force'])
today = datetime.utcnow().strftime("%y%m%d%H%M")
print('Date: {}'.format(today))
Output = "{0}_{1}".format(rPar['output'],today)
mDDC10.runAcq(nFiles=rPar['nRuns'],outDir=Output)
os.system('cp {0} {1}/{2}/{2}_conf.yaml'.format(args.conf,rPar['datadir'],Output))
mDDC10.tn.close()
