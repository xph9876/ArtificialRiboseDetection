#!/usr/bin/env python3

import argparse
import sys
from checkRes import check_res
from subtractionUtils import analyze
import numpy as np
import xlsxwriter

# create worksheets for output
def create_worksheets(filename, re_all, specials):
    workbook = xlsxwriter.Workbook(filename, {'nan_inf_to_errors':True})
    percent = workbook.add_format({'num_format':'0.00%'})
    chroms = ['All', 'Nuc'] + specials
    worksheets = {x:workbook.add_worksheet(x) for x in chroms}
    # output header
    for k, w in worksheets.items():
        w.write(0, 0, 'Sample')
        w.write(0, 1, f'Total_{k}')
        w.write(0, 2, f'No-dA_{k}')
        w.write(0, 3, f'dA_{k}')
        i = 4
        for e in list(re_all.values())[0]:
            w.write(0, i, f'{e}_No-dA')
            w.write(0, i+1, f'{e}_dA')
            i += 2
    return workbook, worksheets, percent


# write results into workbook
def write_workbook(sample, out, worksheets, re_all, row, formatter):
    res = ['all'] + list(re_all.values())[0]
    for k, v in out.items():
        w = worksheets[k]
        w.write(row, 0, sample)
        w.write(row, 1, sum(v['all']), formatter)
        for e in res:
            idx = res.index(e) * 2 + 2
            if e not in v:
                v1 = (np.nan, np.nan)
            else:
                v1 = v[e]
            w.write(row, idx, v1[0], formatter)
            w.write(row, idx+1, v1[1], formatter)


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Calculate ribose ratio in middle positions of restriction enzyme cut sites.')
    parser.add_argument('list', type=argparse.FileType('r'), default=sys.stdin, help='list of calculation, fs\tspecies\tRE1,RE2,....')
    parser.add_argument('raw_reads', help='folder of raw reads')
    parser.add_argument('-o', type=argparse.FileType('wb'), default='RE_calc_results.xlsx', help='output XLSX filename')
    parser.add_argument('-d', default=11, type=int, help='Barcode plus umi length (11)')
    parser.add_argument('--subtracted_output', default=None, help='Folder to output subtracted files. Do not subtract by default')
    parser.add_argument('--bed', default='.', help='Folder of bed files, default=current folder')
    parser.add_argument('--special', type=str, nargs='*', default=['chrM'], help='Special sequences to be calculated, default=chrM')
    parser.add_argument('--RElist', type=argparse.FileType('r'), default='res_all.list', help='pool of RElist, Name\tPattern\tPos, default=res_all.list')
    parser.add_argument('--genome', default=None, help='Folder of genome, provide if needed')
    args = parser.parse_args()

    if args.subtracted_output == '':
        args.subtracted_output = args.bed

    # get res information for libs and generate file
    # libinfo[lib] = [species,[RE1,RE2,...]]
    # re_all[sp] = [REs]
    # residue[re] = ['TNNN']
    libinfo, re_all, residue, da, noda = check_res(args.list, args.RElist, args.bed, args.genome)
    print(da, noda)

    print('Data read!')
    print('Libraries\tREs')
    for k in sorted(libinfo.keys()):
        print(k + '\t'+','.join(libinfo[k][1]))
    print('Bed folder : {}'.format(args.bed))
    print('Raw reads folder : {}'.format(args.raw_reads))
    print('Special sequences : {}'.format(', '.join(args.special)))

    # create output file
    wb, worksheets, formatter = create_worksheets(args.o, re_all, args.special)
    
    # analyze
    # analyze(fs, libinfo, residue, da, noda, bed_folder, fq_folder, d, specials, subtract):
    row = 1
    for lib in libinfo.keys():
        out = analyze(lib, libinfo, residue, da, noda, args.bed, \
            args.raw_reads, args.d, args.special, args.subtracted_output)
        # output to workbook
        write_workbook(lib, out, worksheets, re_all, row, formatter)
        row += 1
        print(f'{lib} Done!')
    wb.close()

    print('Done!')


if __name__ == '__main__':
    main()
