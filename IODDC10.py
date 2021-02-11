import getpass
import telnetlib
import time
import yaml

import subprocess, os
from datetime import date

class IODDC10:
	def __init__(self,SMBIP="192.168.1.50",DDCIP="192.168.1.3",nSam=8192,nEvs=10000,chMask='0x2',dataDir='/data/share/',password=None):
		self.SMBIP = SMBIP
		self.DDCIP = DDCIP
		self.RFA = False
		self.nSam = nSam
		self.nEvs = nEvs
		self.chMask = chMask
		self.dataDir = dataDir
		self.tn = telnetlib.Telnet(self.DDCIP)
		print(self.tn.read_until(b'commands.',timeout=3).decode('ascii'))
		if password is None:
			self.password = getpass.getpass()
		else:
			self.password = password

	@classmethod
	def from_yml(cls,path):
		rPar = {}
		with open(path) as file:
			rPar = yaml.load(file,Loader=yaml.FullLoader)
		return cls(SMBIP=rPar['smbip'],DDCIP=rPar['ddcip'],nSam=rPar['nsam'],nEvs=rPar['ntrg'],chMask=hex(int(rPar['chmask'],2)),dataDir=rPar['datadir'],password=rPar['passwd'])

	def setupDDC10(self,fade=6,force=False):
		print("NOTE:::Restarting samba NEEDS to run as root")
		self.tn.write(b"ls /mnt/share\n")
		pout = self.tn.read_eager().decode('ascii')
		pout += self.tn.read_until(b'root:/>',timeout=3).decode('ascii')
		if "No such file or directory" in pout or force:
			if force:
				os.system('service smb restart')
			self.tn.write(b"mkdir -p /mnt/share\n")
			print("directory made")
			if self.password:
				print("mounting samba")
				thiscmd = "smbmount //{}/SHARE /mnt/share -o username=lzer,password={}\n".format(self.SMBIP,self.password)
				self.tn.write(thiscmd.encode('ascii'))
				passtst = self.tn.read_until(b"----",timeout=3).decode('ascii')
				if "failed" not in passtst:
					self.RFA = True
					print("Ready for Acquisition")
				else:
					print("SMB Failed somehow")
					print(passtst)
			else:
				print("Invalid password Set")

			self.tn.write(b"ls /mnt/share\n")
		else:
			print("smb mount dir already created")
			self.RFA = True
		cmd = "fadc_AD9648 {}".format(fade)
		self.tn.write(cmd.encode('ascii')+b'\n')
		#self.tn.write(b"echo -------\n")

		print(self.tn.read_until(b"----",timeout=3).decode('ascii'))
		self.tn.read_eager()

	def runAcq(self,nFiles=1,outDir='data',delay=0):
		if self.RFA:
			self.tn.write("mkdir -p /mnt/share/{}\n".format(outDir).encode('ascii'))
			self.tn.write("ls /mnt/share/{}\n".format(outDir).encode('ascii'))
			print(self.tn.read_until(b"----",timeout=3).decode('ascii'))
			print(self.tn.read_eager().decode('ascii'))
			for i in range(nFiles):
				self.runAcq("{0}/{0}_{1}".format(outDir,i))
				print("Completed file {}".format(i))
				time.sleep(delay)
		else:
			print("Not Ready For Acquisition!!!!\nHAVE YOU RUN SETUP?\n")
