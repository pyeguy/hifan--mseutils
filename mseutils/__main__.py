from .utils import *
from . import file_parsers
import argparse
import os
from collections import namedtuple, defaultdict
import time

rep_file_tup = namedtuple("rep_file_tup",('rep_fname', 'frag_fname', 'wellid','path'))

class FileWriter(object):
    def __init__(self,fname,headers=None,delim='\t'):
        self.fname = fname
        self.headers = headers
        self.delim = delim
        self.fout = None

    def _open_file(self):
        self.fout = open(self.fname,'w')
        self._write_headers()

    def _write_headers(self):
        if self.headers is not None:
            for header in headers:
                print("{}{}".format(header,self.delim),file=self.fout,end='',sep='')
            print('',file=self.fout)
    
    def __enter__(self):
        self._open_file()
    
    def __exit__(self):
        self.fout.close()

    def close(self):
        self.fout.close()
    
    def write(self,d):
        if self.headers is None:
           self.headers = list(d.keys())
        if not self.fout:
            self._open_file()

        for header in self.headers:
            try:
                pritn("{}{}".format(d[header],self.delim),file=self.fout,end='',sep='')
            except KeyError:
                pritn("{}{}".format('NaN',self.delim),file=self.fout,end='',sep='')
        print('',file=self.fout)




def parse_args():
    parser = argparse.ArgumentParser(prog="MseUtils")
    parser.add_argument("-s","--source_frags", help='combine specs on source_frags, should be a path to folder of replicated data')
    args = parser.parse_args()
    return args


def gen_rep_file_pairs(path):
    stuff = os.scandir(path)
    iddict = defaultdict(dict)
    csv_files = [fname for fname in stuff if fname.lower().endswith(".csv")]
    for fname in csv_files:
        lfname = fname.lower()
        if lfname.endswith('replicated.csv'): 
            sampid = lfname.strip("_replicated.csv")
            iddict[sampid]['rep']=fname
        elif lfname.endswith('rep_fragments.csv'): 
            sampid = lfname.strip("_rep_fragments.csv")
            iddict[sampid]['frag']=fname

    rep_file_pairs = []
    for wellid,fnamed in iddict.items():
        if fnamed['rep'] and fnamed['frag']:
            rft = rep_file_tup(
                os.path.join(path,fnamed['rep']),
                os.path.join(path,fnamed['frag']),
                wellid,
                path,
                )
            rep_file_pairs.append(rft)       
    return rep_file_pairs


def load_and_src_frag(rft):
    mses = file_parsers.load_rep_and_frags_csv(rft.rep_fname,rft.frag_fname)
    combined = src_frags(mses)
    filewriter = FileWriter(os.path.join(rft.path,"{wellid}_combined.{ext}".format(wellid=rft.wellid,ext='.txt')))
    for comb_mse in combined:
        filewriter.write(comb_mse.to_dict())



def main():
    st = time.time()
    args = parse_args()
    if args.source_frags:
        rfts = gen_rep_file_pairs(args.source_frags)
        for rft in tqdm(rfts,desc='Combining Specs'):
            load_and_src_frag(rft)
    et = time.time() - st
    print('Done!')
    pritn('Took {:.2f}min'.format(et/60))




if __name__ == '__main__':
    main()