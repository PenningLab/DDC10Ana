import getpass
import telnetlib
import time
from datetime import date

class IODDC10:
	def __init__(self,HOST="192.168.1.3",nSam=8192,nEvs=10000,chMask='0x1'):
		self.HOST = HOST
		self.RFA = False
		self.nSam = nSam
		self.nEvs = nEvs
		self.chMask = chMask
		self.tn = telnetlib.Telnet(self.HOST)

	def setupDDC10(self,fade=5):
		self.tn.write(b"mkdir -p /mnt/share\n")

		self.password = getpass.getpass()

		if self.password:
			self.tn.write(b"smbmount //192.168.1.2/SHARE /mnt/share -o username=lzer\n")
			self.tn.read_until(b"Password: ")
			self.tn.write(self.password.encode('ascii') + b"\n")
			passtst = self.tn.read_until(b"SMB connection failed").decode('ascii')
			if len(passtst)==0:
				self.RFA = True
				print("Ready for Acquisition")
			else:
				print(passtst)
		else:
			print("Invalid password Set")

		self.tn.write(b"ls /mnt/share\n")
		self.tn.write("fadc_AD9648 {}\n".format(fade).encode('ascii'))

		print(self.tn.read_lazy().decode('ascii'))

	def runAcq(self,outFile='data'):
		if self.RFA:
			cmd = "time /mnt/share/binaries/DDC10_BinWaveCap_ChSel {0} {1} {2} /mnt/share/{3}.bin >> /mnt/share/{3}.log\n".format(self.chMask,self.nSam,self.nEvs,outFile)
			runStart = time.time()
			self.tn.write(cmd.encode('ascii'))
			self.tn.read_until(b"sys")
			with open(outFile+".log",'a') as logfile:
				logfile.write(str(runStart)+"\n")
			print(self.tn.read_lazy().decode('ascii'))
		else:
			print("Not Ready For Acquisition!!!!\nHAVE YOU RUN SETUP?\n")

	def loopAcq(self,nFiles=5,outDir='data'):
		if self.RFA:
			for i in range(nFiles):
				self.runAcq("{0}/{1}\n".format(outDir,i))
				print("Completed file {}".format(i))
		else:
			print("Not Ready For Acquisition!!!!\nHAVE YOU RUN SETUP?\n")

def DoExp(nSam,nEvs,chMask,nFiles,OutName):
	mDDC10 = IODDC10(nSam=nSam,nEvs=nEvs,chMask=chMask)
	mDDC10.setupDDC10()
	today = date.today()
	mDDC10.loopAcq(nFiles=nFiles,outDir="{0}_{1}_{2}_samples_{3}_events".format(OutName,today.strftime("%y%m%d"),nSam,nEvs))
	self.tn.close()
