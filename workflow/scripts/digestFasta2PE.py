#!/usr/bin/env python
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fa', default = '/dev/stdin', help = 'fasta input (default stdin)')
    parser.add_argument('--cutsite', default = 'CATG', help = 'cut site for digestion enzyme (default "CATG")')
    parser.add_argument('--outpref', default = 'digestFasta2PE', help = 'output file name prefix')
    parser.add_argument('--minlen', default = 200, type = int, help = 'minimum fragment length to output (default 200)')
    parser.add_argument('--maxFragUse', default = 10, type = int, help = 'max number of times to use a digest fragment (default 10)') 
    args = parser.parse_args()
    workingseq = []
    cutsite = args.cutsite.upper()
    o1 = open(args.outpref + '_R1.fa','w')
    o2 = open(args.outpref + '_R2.fa','w')
    for line in open(args.fa):
        if line[0] == '>':
            if workingseq != []:
                fraguse = defaultdict(lambda: 0)
                workingseq = "".join(workingseq)
                frags = workingseq.split(cutsite)
                for i,frag1 in enumerate(frags):
                    for k,frag2 in enumerate(frags):
                        if i != k and len(frag1) > args.minlen and len(frag2) > args.minlen and k > i and fraguse[i] <= args.maxFragUse and fraguse[k] <= args.maxFragUse:
                            sn = f'{seq}:{i}:{k}'
                            o1.write('>' + sn + ' 1\n' + frag1 + '\n')
                            o2.write('>' + sn + ' 2\n' + frag2 + '\n')
                            fraguse[i] += 1
                            fraguse[k] += 1
                        if fraguse[i] >= args.maxFragUse:
                            break
            seq = line.strip()[1:].split()[0]
            workingseq = []
        else:
            workingseq.append(line.strip().upper())
    o1.close()
    o2.close()

if __name__ == "__main__": main()