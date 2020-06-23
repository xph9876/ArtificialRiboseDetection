#!/usr/bin/env python3

import argparse
import sys
import queue
import threading
import time
from collections import defaultdict
from checkRes import check_res
from subtractionThreads import SubtractionThread
import numpy as np

# set global thread parameter
exitFlag = False
qLock = threading.Lock()

# argparse
parser = argparse.ArgumentParser(description='Calculate ribose ratio in middle positions of restriction enzyme cut sites.')
parser.add_argument('list', type=argparse.FileType('r'), default=sys.stdin, help='list of calculation, fs\tspecies\tRE1,RE2,....')
parser.add_argument('raw_reads', help='folder of raw reads')
parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='output file')
parser.add_argument('--subtracted_output', default='', help='Folder to output subtracted files, default=same to bed')
parser.add_argument('--ns', action='store_true', help='Do not do subtraction')
parser.add_argument('--bed', default='.', help='Folder of bed files, default=current folder')
parser.add_argument('--special', type=str, nargs='+', default=['chrM'], help='Special sequences to be calculated, default=chrM')
parser.add_argument('--RElist', type=argparse.FileType('r'), default='res_all.list', help='pool of RElist, Name\tPattern\tPos, default=res_all.list')
parser.add_argument('--genome', default=None, help='Folder of genome, provide if needed')
parser.add_argument('--threads', type=int, default=1, help='Number of threads')
args = parser.parse_args()

if args.subtracted_output == '':
    args.subtracted_output = args.bed

# get res information for libs and generate file
# libinfo[lib] = [species,[RE1,RE2,...]]
# re_all[sp] = [REs]
# residue_da[re] = ['TNNN']
lib_info, re_all, residue_da = check_res(args.list, args.RElist, args.bed, args.genome)
print('Data read!')
print('Libraries\tREs')
for k in sorted(lib_info.keys()):
    print(k + '\t'+','.join(lib_info[k][1]))
print('Bed folder : {}'.format(args.bed))
print('Raw reads folder : {}'.format(args.raw_reads))
print('Special sequences : {}'.format(', '.join(args.special)))

# set data
# data[lib][noda/da][re] = 1
data = defaultdict( lambda : defaultdict( lambda : defaultdict(int) ) )

# set queue
workQueue = queue.Queue()
job_types=['total','noda','da']
if not args.ns:
    job_types.append('subtract')
for lib,v in lib_info.items():
    qLock.acquire()
    for s in job_types:
        workQueue.put((lib, v[0], v[1], s))
    qLock.release()
print('Job queue set!')

# set thread
threads = []
for idx in range(args.threads):
    thread = SubtractionThread(idx, workQueue, data, residue_da, \
            args.raw_reads, args.bed, args.subtracted_output, args.special, qLock)
    thread.start()
    threads.append(thread)
print('All threads set! Total threads: {:d}'.format(len(threads)))

# check queue
processed = [False] * len(threads)
while not all(processed):
    for t in range(args.threads):
        if not processed[t]:
            if not threads[t].isAlive():
                if threads[t].exitcode == 1:
                    sys.exit(threads[t].exc_traceback)
                else:
                    processed[t] = True
    time.sleep(1)
print('Calculation finished!')

# output
# header
for k, v in re_all.items():
     # output header
    args.o.write('Library\tSpecies\tRE\tnoda_all\tda_all')
    for i in v:
        args.o.write('\t{}_noda\t{}_da'.format(i,i))
    for sseq in args.special:
        args.o.write(f'\tLibrary\tSpecies\tRE\tnoda_{sseq}\tda_{sseq}')
        for i in v:
            args.o.write('\t{}_noda\t{}_da'.format(i,i))
    args.o.write('\tLibrary\tSpecies\tRE\tnoda_nuc\tda_nuc')
    for i in v:
        args.o.write('\t{}_noda\t{}_da'.format(i,i))
    args.o.write('\n')

# data
for fs in lib_info.keys():
    if lib_info[fs][0] == k:
        # total
        da_chr = data[fs]['da']
        noda_chr = data[fs]['noda']
        total_chr = data[fs]['total']['total']
        if total_chr == 0:
            total_chr = np.nan
        args.o.write('{}\t{}\t'.format(fs, k) + ','.join(lib_info[fs][1]) + '\t{:f}\t{:f}'.format(noda_chr['total']/total_chr,da_chr['total']/total_chr))
        for enzyme in v:
            args.o.write('\t{:f}\t{:f}'.format(noda_chr[enzyme]/total_chr, da_chr[enzyme]/total_chr))
        args.o.write('\t')
        # specs
        for sseq in args.special:
            da_chr = data[fs]['da_'+sseq]
            noda_chr = data[fs]['noda_'+sseq]
            total_chr = data[fs]['total']['total_'+sseq]
            if total_chr == 0:
                total_chr = np.nan
            args.o.write('{}\t{}\t'.format(fs, k) + ','.join(lib_info[fs][1]) + '\t{:f}\t{:f}'.format(noda_chr['total']/total_chr,da_chr['total']/total_chr))
            for enzyme in v:
                args.o.write('\t{:f}\t{:f}'.format(noda_chr[enzyme]/total_chr, da_chr[enzyme]/total_chr))
            args.o.write('\t')
        # nuc
        da_chr = data[fs]['da_n']
        noda_chr = data[fs]['noda_n']
        total_chr = data[fs]['total']['total_n']
        if total_chr == 0:
            total_chr = np.nan
        args.o.write('{}\t{}\t'.format(fs, k) + ','.join(lib_info[fs][1]) + '\t{:f}\t{:f}'.format(noda_chr['total']/total_chr,da_chr['total']/total_chr))
        for enzyme in v:
            args.o.write('\t{:f}\t{:f}'.format(noda_chr[enzyme]/total_chr, da_chr[enzyme]/total_chr))
        args.o.write('\t')
        args.o.write('\n')

print('Done!')


