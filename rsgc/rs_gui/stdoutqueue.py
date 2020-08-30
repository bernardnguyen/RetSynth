import sys, multiprocessing
from multiprocessing.queues import Queue

class StdoutQueue(Queue):
	def __init__(self, log_level='RetSynth:\t'):
		super(StdoutQueue, self).__init__(ctx=multiprocessing.get_context())
		self.logfile = None
		self.log_level = log_level

	def write(self,msg):
		self.put(msg)
		if self.logfile and msg:
			msg = msg.replace('\n','')
			with open(self.logfile, 'a') as log:
				log.write(self.log_level + msg + "\n")

	def flush(self):
		sys.__stdout__.flush()

	def setLogFile(self, file):
		self.logfile = file
