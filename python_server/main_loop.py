import tornado.ioloop
import tornado.web
from tornado.log import enable_pretty_logging
import svgwrite
import sys
from falcon_kit import kup, falcon, DWA, get_alignment
from falcon_kit.FastaReader import FastaReader
import subprocess
import shlex
import logging
import numpy as np
import matplotlib as plt
from pylab import *

enable_pretty_logging()

rmap = dict(zip("ACGTN","TGCAN"))


def get_kmer_matches(seq1, seq0):
    K = 8
    lk_ptr = kup.allocate_kmer_lookup( 1 << (K * 2) )
    sa_ptr = kup.allocate_seq( len(seq0) )
    sda_ptr = kup.allocate_seq_addr( len(seq0) )
    kup.add_sequence( 0, K, seq0, len(seq0), sda_ptr, sa_ptr, lk_ptr)
    #kup.mask_k_mer(1 << (K * 2), lk_ptr, 16)
    kmer_match_ptr = kup.find_kmer_pos_for_seq(seq1, len(seq1), K, sda_ptr, lk_ptr)
    kmer_match = kmer_match_ptr[0]
    aln_range = kup.find_best_aln_range(kmer_match_ptr, K, K*50, 50)
    
    x,y = zip( * [ (kmer_match.query_pos[i], kmer_match.target_pos[i]) for i in range(kmer_match.count )] )
    kup.free_kmer_match(kmer_match_ptr)
    return x, y, aln_range



class Mode1(tornado.web.RequestHandler):

    def get(self, in_data):
    	print in_data
    	self.write("Hello, world")

    def post(self):
        chr_ = self.get_argument("chr", "NA")
        start = int(self.get_argument("start", -1))
        end = int(self.get_argument("end", -1))
        strand = self.get_argument("strand","NA")
        read_name = self.get_argument("read_name", "NA")
        read_seq = self.get_argument("read_seq", "NA")
        ref_seq = self.get_argument("ref_seq", "NA")

        

        #the sequnce from IGV is always the same direction as the reference  
        #if strand == "NEGATIVE":
        #	read_seq = "".join( [rmap[c] for c in read_seq[::-1]] )

    	x, y, ar = get_kmer_matches(read_seq, ref_seq)
        y = np.array(y) + start
        plt.figure()
        plt.plot(x, y, "o", markeredgecolor="none", alpha=0.5)
        xlabel( read_name )
        ylabel( chr_ )

        plt.show()
    	self.write("Dot Plot Done") 

class MainHandler(tornado.web.RequestHandler):
	def get(self):
		self.write("Hello, world")

application = tornado.web.Application([
    (r"/", MainHandler),
    (r"/mode1/", Mode1),], 
    autoreload=True, debug=True)

if __name__ == "__main__":
    from tornado.options import options
    options.logging = "DEBUG"
    logging.debug("starting torando web server")
    application.listen(6502)
    tornado.ioloop.IOLoop.current().start()
