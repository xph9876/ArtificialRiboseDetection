import sys
import subprocess
import re
import os
from collections import defaultdict
from reBuild import re_build

# A part of calculate percentage of ribos incorporated in restriction enzyme cut site
# TODO: 1. get all restriction enzyme used by libaries
      # 2. generate cuting pattern with dA-tailing case
      # 3. check whether bed folder contains bed file of res cutting site
# INPUT: lib: filereader for libraries
#        re_list: filereader for RE list
#        folder : folder path for bed files
#        genome_folder : folder path for genome files
def check_res(lib, re_list, folder, genome_folder):
    # build re information for libaries used
    # libinfo[lib]=[species,[res]]
    # re_all[species] = [res]
    libinfo = {}
    re_all = defaultdict(list)
    for l in lib:
        ws = l.rstrip('\n').split('\t')
        if len(ws) != 3:
            continue
        res = [x.lower() for  x in ws[-1].split(',')]
        libinfo[ws[0]] = [ws[1], res]
        # add re to all list
        for enzyme in res:
            if enzyme not in re_all[ws[1]]:
                re_all[ws[1]].append(enzyme)

    # find cut pattern
    degenerated_base = {'W':'[AT]', 'S':'[CG]', 'M':'[AC]', 'K':'[GT]', 'R':'[AG]', 'Y':'[CT]',\
            'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'}
    residue_da = {}
    pattern = {}
    allre = []
    for k,v in re_all.items():
        allre += v
    for l in re_list:
        ws = l.split('\t')
        ws[0] = ws[0].lower()
        if len(ws) !=3:
            continue
        if ws[0] in allre:
            # calc_residue
            residue = 'T' + ws[1][int(ws[2]):]
            # replace
            for k,v in degenerated_base.items():
                ws[1] = ws[1].replace(k,v)
                residue = residue.replace(k,v)
            pattern[ws[0]] = ws[1]
            residue_da[ws[0]] = residue

    # check whether exist da and noda
    for sp, v in re_all.items():
        for enzyme in v:
            if not (os.path.exists(folder + '/' + '_'.join([sp, enzyme, 'noda.bed'])) and \
                    os.path.exists(folder + '/' + '_'.join([sp, enzyme, 'da.bed']))):
                print('Cannot find da/noda bed file for {}_{} in folder {}, generate it!'.format(sp, enzyme, folder))
                if not genome_folder:
                    sys.exit('[ERROR]Cannot find genome folder, please provide it in parameters')
                try:
                    with open(genome_folder + '/{}.fa'.format(sp)) as genome:
                        re_build(pattern[enzyme], genome, folder + '/{}_{}'.format(sp,enzyme))
                except IOError:
                    sys.exit('[ERROR]Cannot open genome {}.fa at folder{}!'.format(sp, genome_folder))

    return(libinfo, re_all, residue_da)


