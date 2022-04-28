import regex as re


# find lenth of pattern
def determine_length(pattern):
    lenc = 0
    i = 0
    skip = False
    while i< len(pattern):
        if pattern[i] == '[':
            lenc += 1
            skip = True
            i += 1
        elif pattern[i] == ']':
            skip = False
            i += 1
        elif pattern[i] == '{':
            i += 1
            num = ''
            while pattern[i] != '}':
                num += pattern[i]
                i += 1
            lenc += (int(num)-1)
            i += 1
        elif skip:
            i += 1
        else:
            lenc += 1
            i += 1
    lenc -= 1
    return lenc
 

def re_build(patterns, enzymes, fr):
    # find max pattern length as cache
    lengths = {x:determine_length(patterns[x]) for x in enzymes}
    lenc = max(lengths.values())

    # find da site and noda site
    da = {x:{} for x in enzymes}
    noda = {x:{} for x in enzymes}
    for l in fr:
        l = l.rstrip('\n')
        if l[0] == '>':
            chrom = l[1:]
            cache = ''
            pos = 0
            assert chrom not in da, f'Duplicate chromosome {chrom} found in fasta file'
            for e in enzymes:
                da[e][chrom] = set()
                noda[e][chrom] = set()
        else:
            l = cache + l.upper()
            for e in enzymes:
                pattern = patterns[e]
                matches = re.finditer(pattern,l, overlapped=True)
                for m in matches:
                    mid = pos + (m.start() + m.end())//2
                    noda[e][chrom].add((mid-1, '+'))
                    noda[e][chrom].add((mid, '-'))
                    da[e][chrom].add((mid-1, '-'))
                    da[e][chrom].add((mid, '+'))
            cache = l[-lenc:]
            pos = pos + len(l) -lenc
    
    return noda, da
