import re

def re_build(pattern, fr, base):
    # find length of find pattern
    lenc = 0
    i = 0
    skip = False
    while i< len(pattern):
        if pattern[i] == '[':
            lenc += 1
            skip=True
            i += 1
            continue
        elif pattern[i] == ']':
            skip = False
            i += 1
            continue
        elif pattern[i] == '{':
            i += 1
            num = ''
            while pattern[i] != '}':
                num += pattern[i]
                i += 1
            lenc += (int(num)-1)
            i += 1
            continue
        elif skip:
            i += 1
            continue
        else:
            lenc += 1
            i += 1
    lenc -= 1

    # find da site and noda site
    with open(base + '_da.bed', 'w') as da, open(base + '_noda.bed', 'w') as noda:
        for l in fr:
            l = l.rstrip('\n')
            if l== '':
                continue
            elif l[0] == '>':
                chrom = l[1:]
                cache = ''
                pos = 0
            else:
                l = cache + l.upper()
                matches = re.finditer(pattern,l)
                for m in matches:
                    mid = int(pos + (m.start() + m.end())/2)
                    noda.write('\t'.join([chrom, str(mid-1), str(mid), '.', '.', '+']) + '\n')
                    noda.write('\t'.join([chrom, str(mid), str(mid+1), '.', '.', '-']) + '\n')
                    da.write('\t'.join([chrom, str(mid), str(mid+1), '.', '.', '+']) + '\n')
                    da.write('\t'.join([chrom, str(mid-1), str(mid), '.', '.', '-']) + '\n')
                cache = l[-lenc:]
                pos = pos + len(l) -lenc
