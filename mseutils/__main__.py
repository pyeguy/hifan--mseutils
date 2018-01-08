from .utils import *
from . import file_parsers
from . import mseh5
from . import msesql
import argparse
import os
from collections import namedtuple, defaultdict
import time
from concurrent.futures import ProcessPoolExecutor,as_completed


rep_file_tup = namedtuple("rep_file_tup",('rep_fname', 'frag_fname', 'sampid','path'))

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
        for header in self.headers:
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
                print("{}{}".format(d[header],self.delim),file=self.fout,end='',sep='')
            except KeyError:
                print("{}{}".format('NaN',self.delim),file=self.fout,end='',sep='')
        print('',file=self.fout)

def parse_args():
    parser = argparse.ArgumentParser(prog="MseUtils")
    parser.add_argument("-s","--source_frags", help='combine specs on source_frags, should be a path to folder of replicated data')
    parser.add_argument("-p","--procs",help='number of procs to use on multiproc functions',default=2,type=int)
    parser.add_argument("--tolerances",help='print the tolerances used',action="store_true")
    parser.add_argument("--h5",help="write an HDF5 file with all the Mse data in the path",action="store_true")
    parser.add_argument("--sqlite",help="write to a SQLite database with all the Mse data",action="store_true")
    
    args = parser.parse_args()
    return args

def gen_rep_file_pairs(path):
    stuff = os.scandir(path)
    iddict = defaultdict(dict)
    files = [os.path.basename(name) for name in stuff if os.path.isfile(name)]
    csv_files = [fname for fname in files if fname.lower().endswith(".csv")]
    for fname in csv_files:
        lfname = fname.lower()
        if lfname.endswith('replicated.csv'): 
            sampid = fname.strip("_replicated.csv")
            iddict[sampid]['rep']=fname
        elif lfname.endswith('rep_fragments.csv'): 
            sampid = fname.strip("_rep_fragments.csv")
            iddict[sampid]['frag']=fname

    rep_file_pairs = []
    for sampid,fnamed in iddict.items():
        try:
            rft = rep_file_tup(
                os.path.join(path,fnamed['rep']),
                os.path.join(path,fnamed['frag']),
                sampid,
                path,
                )
            rep_file_pairs.append(rft)       
        except KeyError:
            # add verbose conditional
            print("Bad sampid {}".format(sampid))
    return rep_file_pairs

def load_and_src_frag(rft,write_flat_file=False):
    mses = file_parsers.load_rep_and_frags_csv(rft.rep_fname,rft.frag_fname)
    combined = src_frags(mses)
    sampid = rft.sampid
    if write_flat_file:
        filewriter = FileWriter(os.path.join(rft.path,"{sampid}_combined.{ext}".format(sampid=rft.sampid,ext='txt')))
        for comb_mse in combined:
            for fragd in comb_mse.to_frags():
                filewriter.write(fragd)
        return None
    else:
        return sampid,combined

def clearscreen(numlines=100):
    """Clear the console.
    numlines is an optional argument used only as a fall-back.
    """
    # Thanks to Steven D'Aprano, http://www.velocityreviews.com/forums

    if os.name == "posix":
        # Unix/Linux/MacOS/BSD/etc
        os.system('clear')
    elif os.name in ("nt", "dos", "ce"):
        # DOS/Windows
        os.system('CLS')
    else:
    # Fallback for other operating systems.
        print('\n' * numlines)

def main():
    st = time.time()
    args = parse_args()
    if args.tolerances:
        print_tolerances()
    if args.source_frags:
        rfts = gen_rep_file_pairs(args.source_frags)
        if args.h5:
            h5t = mseh5.create_h5_file("Output.h5")
        if args.sqlite:
            conn = msesql.create_db("Output")

        ################################
        ## Single Proc for Debugging ##
        ################################
        # idx =0
        # for rft in rfts:
        #     sampid,mses = load_and_src_frag(rft)
        #     if args.h5:
        #          idx = mseh5.add_mses(mses,h5t,sampid=sampid,idx=idx)
 
        with ProcessPoolExecutor(max_workers=args.procs)  as executor:
            # clearscreen()
            pbar = tqdm(total=len(rfts),desc='Combining Specs')
            if not args.h5 and not args.sqlite:
                futs = [executor.submit(load_and_src_frag,rft,write_flat_file=True) for rft in rfts]
            else:
                futs = [executor.submit(load_and_src_frag,rft) for rft in rfts]
            idx =0
            for fut in as_completed(futs):
                if args.h5 or args.sqlite:
                    e = fut.exception()
                    if e:
                        raise e
                    sampid,mses = fut.result()
                    if args.h5:
                        idx = mseh5.add_mses(mses,h5t,sampid=sampid,idx=idx)
                    if args.sqlite:
                        idx = msesql.add_mses(conn,mses,sampid=sampid,idx=idx)

                pbar.update()
        pbar.close()
        if args.h5:
            h5t.h5.close()
    et = time.time() - st
    print('Done!')
    print('Took {:.2f}min'.format(et/60))

if __name__ == '__main__':
    main()