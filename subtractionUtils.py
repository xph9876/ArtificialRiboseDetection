import subprocess
import os
import re

# check read if started with da
def check_read(s, pattern):
    # ribose start from 11th pos
    if re.match(pattern, s[11:]):
        return True
    else:
        return False


# get reads from a list and raw fastq file
def find_reads_name(nd, pattern, fr, names):
    captured = False
    for l in fr:
        if captured:
            # check if a da-tailing self ligation
            if check_read(l, pattern):
                names[name + nd[name]] = 1
            captured = False
        # search name
        elif l[0] == '@':
            name = l[1:].split(' ')[0]
            if name in nd:
                captured = True


# get reads from a dict and raw fastq file
def count_da(names, pattern, fr):
    da_ribo = 0
    captured = False
    for l in fr:
        if captured:
            # check if a da-tailing self ligation
            if check_read(l, pattern):
                da_ribo += 1
            captured = False
        # search name
        elif l[0] == '@':
            name = l[1:].split(' ')[0]
            if name in names:
                captured = True
    return da_ribo


# total
def calc_total(data, fs, folder, mitochondrial_name):
    # get total number
    lines = subprocess.Popen(['grep', mitochondrial_name, folder+'/'+fs + '.bed'], stdout=subprocess.PIPE)
    total_m = int(subprocess.run(['wc', '-l'], stdin=lines.stdout, stdout=subprocess.PIPE).stdout.decode(encoding='utf-8'))
    total = int(subprocess.run(['wc', '-l', folder+'/'+fs + '.bed'], stdout=subprocess.PIPE).stdout.decode(encoding='utf-8').split(' ')[0])
    total_n = total - total_m
    data[fs]['total']['total'] = total
    data[fs]['total']['total_n'] = total_n
    data[fs]['total']['total_m'] = total_m


# noda
def calc_noda(data, fs, species, res, folder, mitochondrial_name):
    # noda, use bedtools intersect and then count lines
    noda = data[fs]['noda']
    noda_m = data[fs]['noda_m']
    noda_n = data[fs]['noda_n']
    for enzyme in res:
        inter = subprocess.Popen(['bedtools', 'intersect', '-a', folder+'/'+fs + '.bed', '-b', folder+'/'+'{}_{}_noda.bed'.format(species, enzyme),'-s', '-nonamecheck'],stdout=subprocess.PIPE)
        noda[enzyme] = int(subprocess.run(['wc', '-l'], stdin=inter.stdout, stdout=subprocess.PIPE).stdout.decode(encoding='utf-8'))

        # noda, split mitochondria
        chrm = subprocess.Popen(['grep', mitochondrial_name, folder+'/'+fs + '.bed'], stdout=subprocess.PIPE)
        inter_m = subprocess.Popen(['bedtools', 'intersect', '-a', '-', '-b', folder+'/'+'{}_{}_noda.bed'.format(species, enzyme), '-s', '-nonamecheck'],stdin=chrm.stdout, stdout=subprocess.PIPE)
        noda_m[enzyme] = int(subprocess.run(['wc', '-l'], stdin=inter_m.stdout, stdout=subprocess.PIPE).stdout.decode(encoding='utf-8'))
        noda_n[enzyme] = noda[enzyme] - noda_m[enzyme]
    # total noda
    noda['total'] = sum([ noda[i] for i in res])
    noda_m['total'] = sum([noda_m[i] for i in res])
    noda_n['total'] = sum([noda_n[i] for i in res])


# da
def calc_da(data, fs, species, res, residue_da, rawreads, folder, mitochondrial_name):
    # da
    da = data[fs]['da']
    da_m = data[fs]['da_m']
    da_n = data[fs]['da_n']
    for enzyme in res:

        # get read names for da
        inter = subprocess.Popen(['bedtools', 'intersect', '-a', folder+'/'+fs + '.bed', '-b', folder+'/'+'{}_{}_da.bed'.format(species, enzyme), '-s','-nonamecheck'],stdout=subprocess.PIPE)
        name_l= subprocess.run(['cut', '-f', '4'], stdin=inter.stdout, stdout=subprocess.PIPE).stdout.decode(encoding='utf-8').split('\n')

        names = {}
        # remove adapter information
        for i in range(len(name_l)):
            names[name_l[i].split('_')[0]] = 1

        # check data
        with open(rawreads + '/{}.fq'.format(fs), 'r') as fr:
            da[enzyme] = count_da(names, residue_da[enzyme], fr)

        # da, mitochondria
        chrm = subprocess.Popen(['grep', mitochondrial_name, folder+'/'+ fs + '.bed'], stdout=subprocess.PIPE)
        inter = subprocess.Popen(['bedtools', 'intersect', '-a', '-', '-b', folder+'/'+'{}_{}_da.bed'.format(species, enzyme), '-s','-nonamecheck'],stdin=chrm.stdout, stdout=subprocess.PIPE)
        name_l= subprocess.run(['cut', '-f', '4'], stdin=inter.stdout, stdout=subprocess.PIPE).stdout.decode(encoding='utf-8').split('\n')

        # remove adapter information
        names = {}
        for i in range(len(name_l)):
            names[name_l[i].split('_')[0]] = 1

        # check data
        with open(rawreads + '/{}.fq'.format(fs), 'r') as fr:
            da_m[enzyme] = count_da(names, residue_da[enzyme], fr)
            da_n[enzyme] = da[enzyme] - da_m[enzyme]

    # total da
    da['total'] = sum([da[i] for i in res])
    da_m['total'] = sum([da_m[i] for i in res])
    da_n['total'] = sum([da_n[i] for i in res])


# subtract
def subtract(fs, species, res, residue_da, raw_reads, folder, outfolder, qLock):
    # reads to cut
    name_to_remove = {}
    for enzyme in res:
        # get read names for da
        inter = subprocess.Popen(['bedtools', 'intersect', '-a', folder+'/'+fs+'.bed', '-b', folder+'/{}_{}_da.bed'.format(species, enzyme), '-s','-nonamecheck'],stdout=subprocess.PIPE)
        names= subprocess.run(['cut', '-f', '4'], stdin=inter.stdout, stdout=subprocess.PIPE).stdout.decode(encoding='utf-8').split('\n')

        # remove adapter information
        nd = {}
        for i in names:
            if i == '':
                continue
            a,b = i.split('_')
            nd[a] = '_' + b

        # check data
        with open(raw_reads + '/{}.fq'.format(fs), 'r') as fr:
            find_reads_name(nd, residue_da[enzyme], fr, name_to_remove)

    # cut noda
    comm = ['cat']
    name = folder+'/{}'.format(species)
    for enzyme in res:
        comm += [folder+'/{}_{}_noda.bed'.format(species, enzyme)]
        name += '_{}'.format(enzyme)
    name += '.bed'
    # generate union of noda files
    if not os.path.exists(name):
        qLock.acquire()
        print('Generate union of RE cut site!')
        print('Species: {}\t'.format(species), 'RES: {}'.format(', '.join(res)))
        with open(name,'wb') as nodaw:
            noda = subprocess.Popen(comm, stdout=nodaw)
        qLock.release()
    # subtract noda
    cache_file_name = '.noda_subed_{}.bed'.format(fs)
    with open(cache_file_name, 'w+') as cachew, open(outfolder+'/{}_subtracted.bed'.format(fs), 'w') as fw:
        noda_subed = subprocess.run(['bedtools', 'subtract', '-a', folder+'/'+fs+'.bed', '-b', name, '-s','-nonamecheck'], stdout=cachew)
        cachew.seek(0)
        for l in cachew:
            if not l.split('\t')[3] in name_to_remove:
                fw.write(l)
    # remove cache
    os.remove(cache_file_name)

