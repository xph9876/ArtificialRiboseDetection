import sys
import traceback
import threading
import time
from subtractionUtils import *

# class for threading
class SubtractionThread(threading.Thread):
    def __init__(self, idx, q, data, residue_da, rawreads, bed_folder, outfolder, mitochondrial_name, qLock):
        threading.Thread.__init__(self)
        self.idx = idx
        self.q = q
        self.data = data
        self.folder = bed_folder
        self.outfolder = outfolder
        self.residue_da = residue_da
        self.rawreads = rawreads
        self.mitochondrial_name = mitochondrial_name
        # exception
        self.exitcode = 0
        self.exception = None
        self.exc_traceback = ''
        # function
        self.qLock = qLock
        self.stop_thread = False

    # deal with exception
    def run(self):
        try:
            print('Thread {:d} Start!'.format(self.idx))
            self.get_job(self.data, self.q, self.idx, self.residue_da, self.rawreads, self.folder, self.outfolder, self.mitochondrial_name)
            print('Thread {:d} Closed!'.format(self.idx))
        except Exception as e:
            self.exitcode = 1
            self.exception = e
            self.exc_traceback = '[ERROR] Thread {:d}:\n'.format(self.idx)+''.join(traceback.format_exception(*sys.exc_info()))

    # simple function to choose da or noda
    def get_job(self, data, q, idx, residue_da, rawreads, folder, outfolder, mitochondrial_name):
        # q : (fs, species, res, status)
        while not self.stop_thread:
            self.qLock.acquire()
            if not q.empty():
                # only three status, total, da tailing or no da tailing
                (fs,species, res, status) = q.get()
                self.qLock.release()
                if status == 'total':
                    print ("Thread {:d} is processing {} total!".format(idx, fs))
                    calc_total(data, fs, folder, mitochondrial_name)
                    print ("{} total finished!".format(fs))
                elif status == 'da':
                    print ("Thread {:d} is processing {} da!".format(idx, fs))
                    calc_da(data, fs, species, res, residue_da, rawreads, folder, mitochondrial_name)
                    print ("{} da finished!".format(fs))
                elif status == 'noda':
                    print ("Thread {:d} is processing {} noda!".format(idx, fs))
                    calc_noda(data, fs, species, res, folder, mitochondrial_name)
                    print ("{} noda finished!".format(fs))
                elif status == 'subtract':
                    print('Thread {:d}: Start subtraction for {}!'.format(idx, fs))
                    subtract(fs, species, res, residue_da, rawreads, folder, outfolder, self.qLock)
                    print('Thread {:d}: {} subtraction finished!'.format(idx, fs))
                else:
                    sys.exit('Thread {:d} cannot deal with job type {} of {}!'.format(idx, status, fs))
            else:
                self.stop_thread = True
                self.qLock.release()
            time.sleep(1)

