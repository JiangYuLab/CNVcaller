# -*- coding: utf-8 -*-


import click
from collections import OrderedDict



def loadfa(fafile):
    seqdict = OrderedDict()
    with open(fafile) as f:
        seq_id = f.readline().strip().split()[0][1:]
        tmp_seq = []
        for line in f:
            if line[0] != '>':
                tmp_seq.append(line.strip())
            else:
                seqdict[seq_id] = ''.join(tmp_seq)
                seq_id = line.strip().split()[0][1:]
                tmp_seq = []
        seqdict[seq_id] = ''.join(tmp_seq)
    return seqdict


@click.command()
@click.argument('fafile')
@click.argument('winsize', type=int)
@click.argument('outfile')
def main(fafile, winsize, outfile):
    """
    produce kmer(<OUTFILE>) from refrence fasta(<FAFILE>) according to window size(<WINSIZE>)
    """
    seqdict = loadfa(fafile)
    stepsize = int(winsize / 2)
    print(f'winsize is {winsize}')
    print(f'stepsize is {stepsize}')
    with open(outfile, 'w') as f:
        for name, seq in seqdict.items():
            for start in range(0, len(seq), stepsize):
                f.write(f'>{name}_{start}\n{seq[start: start+winsize]}\n')

if __name__ == '__main__':
    main()

