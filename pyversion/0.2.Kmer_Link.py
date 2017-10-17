# -*- coding: utf-8 -*-


import click
import numpy as np
import pandas as pd


def rename(x):
    """
    reshape seq ID from chr1_0 to chr1-1
    _ -> -
    pos += 1
    """
    id_ = '_'.join(x.split('_')[:-1])
    pos = int(x.split('_')[-1])
    pos += 1
    return f'{id_}-{pos}'



@click.command()
@click.argument('blasr')
@click.argument('winsize', type=int)
@click.argument('outfile')
def main(blasr, winsize, outfile):
    r"""
    produce link file from blasr outputfile, blasr use -m 5 out format
    """
    # blasr outfmt https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format
    stepsize = int(winsize / 2)
    print(f'winsize {winsize}. stepsize {stepsize}')
    with open(blasr) as f:
        for line in f:
            sline = line.strip().split()
            sline[0] = sline[0].split('_')
    col_names = ['qName', 'qLength', 'qStart', 'qEnd', 'tName', 'tStart', 'numMatch', 'numMismatch', 'numIns', 'numDel']
    col_index = [0, 1, 2, 3, 5, 7, 11, 12, 13, 14]
    df = pd.read_csv(blasr, sep='\s+', header=None, usecols=col_index, names=col_names)
    df['qName'] = df['qName'].map(rename)
    print(f'parse blasr finsh, records number is {df.shape[0]}')
    df['coverage'] = (df['qEnd'] - df['qStart']) / df['qLength']
    df = df.loc[df['coverage'] >= 0.9, ['qName', 'qStart', 'tName', 'tStart', 'numMatch', 'numMismatch', 'numIns', 'numDel']]
    print(f'drop coverage lower than 90%, {df.shape[0]} records remained')
    df['ident'] = df['numMatch'] / (df['numMatch'] + df['numMismatch'] + df['numIns'] + df['numDel'])
    df = df.loc[df['ident'] >= 0.97, ['qName', 'tName', 'tStart']]
    print(f'drop identity lower than 97%, {df.shape[0]} records remained')
    df['tNewstart'] = df['tStart'].map(lambda x: round(x / stepsize) * stepsize + 1)
    df['tNewname'] = df['tName'].astype(str) + '-' + df['tNewstart'].astype(str)
    df = df.groupby('qName')['qName', 'tNewname'].agg({'qName': np.unique,
                                                      'tNewname': lambda x: '\t'.join(x)})
    with open(outfile, 'w') as f:
        for line in df.values:
            if len(line[1].split()) > 1:
                outline = [str(x) for x in line]
                f.write('\t'.join(outline) + '\n')


if __name__ == '__main__':
    main()
