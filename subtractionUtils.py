import numpy as np
import re

# check read if started with da pattern
def match_da_pattern(s, pattern, d):
    return re.match(pattern, s[d:]) != None


# update noda with the new line
def update_noda(noda_out, noda, chrom, loc, st, specials, enzymes):
    is_noda = False
    for e in enzymes:
        if (loc, st) in noda[e][chrom]:
            if not is_noda:
                is_noda = True
            # update noda
            noda_out['All'][e] += 1
            if chrom in specials:
                noda_out[chrom][e] += 1
            else:
                noda_out['Nuc'][e] += 1
    # update all enzymes combined
    if is_noda:
        noda_out['All']['all'] += 1
        if chrom in specials:
            noda_out[chrom]['all'] += 1
        else:
            noda_out['Nuc']['all'] += 1
    return is_noda


# update da with the new line
def update_da(da_out, da, chrom, loc, read_name, st, specials, enzymes):
    is_da = False
    for e in enzymes:
        if (loc, st) in da[e][chrom]:
            if not is_da:
                is_da = True
            # update da
            da_out['All'][e].add(read_name)
            if chrom in specials:
                da_out[chrom][e].add(read_name)
            else:
                da_out['Nuc'][e].add(read_name)
    # update all enzymes combined
    if is_da:
        da_out['All']['all'].add(read_name)
        if chrom in specials:
            da_out[chrom]['all'].add(read_name)
        else:
            da_out['Nuc']['all'].add(read_name)
    return is_da
    

# analyze a single library
def analyze(fs, libinfo, residue, da, noda, bed_folder, fq_folder, d, specials, subtract):
    # extract information
    species, enzymes = libinfo[fs]
    # initialization
    total = {x:0 for x in specials + ['All', 'Nuc']}
    noda_out = {x:{y:0 for y in enzymes + ['all']} for x in specials + ['All', 'Nuc']}
    da_out = {x:{y:set() for y in enzymes + ['all']} for x in specials + ['All', 'Nuc']}
    if subtract:
        outbed = []
        da_cache = {}
        fw = open(f'{subtract}/{fs}.bed', 'w')
    # go through bed file to check noda and total lines
    with open(f'{bed_folder}/{fs}.bed', 'r') as fr:
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 6:
                continue
            # get info
            loc = int(ws[1])
            chrom = ws[0]
            st = ws[5]
            read_name = ws[3].split('_')[0]
            # update total
            total['All'] += 1
            if chrom in specials:
                total[chrom] += 1
            else:
                total['Nuc'] += 1
            # update noda
            is_noda = update_noda(noda_out, noda, chrom, loc, st, specials, enzymes)
            # update da
            is_da = update_da(da_out, da, chrom, loc, read_name, st, specials, enzymes)
            # append to output if subtracted files are needed
            if not subtract:
                continue
            if not is_noda and not is_da:
                outbed.append(ws)
            elif not is_noda:
                da_cache[read_name] = ws
    # go through fastq file to confirm da status
    final_da_out = {x:{y:0 for y in enzymes + ['all']} for x in specials + ['All', 'Nuc']}
    with open(f'{fq_folder}/{fs}.fq', 'r') as fr:
        i = -1
        for l in fr:
            i += 1
            # check header
            if not i % 4:
                read_name = l[1:].rstrip('\n').split(' ')[0]
                if read_name not in da_cache:
                    captured = False
                else:
                    captured = True
            # check reads
            elif i % 4 == 1:
                if not captured:
                    continue
                is_da = False
                for chrom, v in da_out.items():
                    is_da_in_chrom = False
                    for e in enzymes:
                        reads = v[e]
                        if read_name not in reads:
                            continue
                        # is dA tailing
                        if match_da_pattern(l, residue[e], d):
                            is_da_in_chrom = True
                            is_da = True
                            final_da_out[chrom][e] += 1
                    if is_da_in_chrom:
                        final_da_out[chrom]['all'] += 1
                if not is_da and subtract != None:
                    outbed.append(da_cache[read_name])
    # output bed
    if subtract:
        outbed.sort()
        for l in outbed:
            fw.write('\t'.join([str(x) for x in l]) + '\n')
        fw.close()
    # calculate final output
    out = {x:{y:0 for y in enzymes + ['all']} for x in specials + ['All', 'Nuc']}
    for chrom, v in noda_out.items():
        for e in v.keys():
            if not total[chrom]:
                out[chrom][e] = (np.nan, np.nan)
            else:
                out[chrom][e] = (noda_out[chrom][e]/total[chrom], final_da_out[chrom][e]/total[chrom])
    return out
