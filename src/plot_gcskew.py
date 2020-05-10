#!/bin/python
#Jennifer Lu, jlu26@jhmi.edu
#2019/01/15
#
#This program calculates and plots gc_skew for a genome provided to this program
import sys, os, argparse 
from time import gmtime
from time import strftime
from Bio import SeqIO
import numpy as np

import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.ticker as tick

def reformat_bp(tick_val, pos):
    val = round(tick_val/1000000,1)
    new_tick = '{:} Mb'.format(val)
    new_tick = str(new_tick)
    return new_tick

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input","-s","--seq", 
        dest="in_file", required=True,
        help="Sequence for which to calculate gc_skew")
    parser.add_argument("-k","-l","--window_len", default=20000,
        dest="window_size", required=False, type=int,
        help="Window size for which to calculate each g-c/g+c [default: 20Kbp]")
    parser.add_argument("--freq",required=False,type=int,
        dest="freq", default=1000, 
        help="Frequency at which to calculate GC skew [default: 1Kbp]")
    parser.add_argument("-o","--output", required=False, default="curr.png", 
        dest="out_file", help="Name of GC Skew Plot to Save to (default: curr.png)")
    args=parser.parse_args() 
    
    #Start program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("   PROGRAM START TIME: " + time + '\n')
    if args.freq == 0:
        freq = args.window_size
    else:
        freq = args.freq
    #Process sequence file 
    id2string = {}
    id2name = {}
    count_seqs = 0 
    sys.stdout.write("\t>> Processing sequence file\n")
    sys.stdout.write("\t\t%i seqs found" % (count_seqs))
    sys.stdout.flush()
    for record in SeqIO.parse(args.in_file,'fasta'):
        count_seqs += 1
        if count_seqs % 100 == 0:
            sys.stdout.write("\r\t\t%i seqs found" % (count_seqs))
            sys.stdout.flush()
        #Save string
        id2name[count_seqs] = record.id
        id2string[count_seqs] = str(record.seq)
    sys.stdout.write("\r\t\t%i seqs found\n" % (count_seqs))
    sys.stdout.flush()

    #Calculate and plot gc skew
    tot = count_seqs
    count_seqs = 0
    sys.stdout.write("\t>> Plotting GC Skew for %i seqs\n" % tot)
    sys.stdout.write("\t\t%i sequences processed" % (count_seqs))
    sys.stdout.flush()
    plot_filename = args.out_file
    fig,ax = plt.subplots(nrows=tot, ncols=1,figsize=(10,3*tot),squeeze=False)
    for i in range(1,tot+1): 
        #Calculate/plot skew
        my_seq = id2string[i]
        my_description = id2name[i]
        count = 0
        #Calculate skew
        g = 0.0
        c = 0.0
        curr_seq = ""
        indices = []
        skews = []
        for j in range(0, len(my_seq), freq):
            if (j+args.window_size) >= len(my_seq):
                break
            curr_seq = my_seq[j:j+args.window_size]
            g = float(curr_seq.count("G"))
            c = float(curr_seq.count("C"))
            if (g+c) > 0:
                new_calc = (g-c)/(g+c)
            else:
                new_calc = 0.0
            #Save values
            indices.append(j)
            skews.append(new_calc)
        #Final print 
        count_seqs += 1
        sys.stdout.write("\r\t\t%i sequences processed" % (count_seqs))
        sys.stdout.flush()
        #Split to pos/neg
        s = np.array(skews)
        pos = np.ma.masked_where(s <= 0, s)
        neg = np.ma.masked_where(s >= 0, s)
        mid = np.ma.masked_where(s != 0, s)
        #Print to file 
        lines = ax[i-1,0].plot(indices,pos, 'deepskyblue', indices,mid, 'g', indices,neg,'k')
        ax[i-1,0].set(xlabel='Genome Position', ylabel='GC Skew',
            title=my_description)
        ax[i-1,0].xaxis.set_major_formatter(tick.FuncFormatter(reformat_bp))
        ax[i-1,0].grid()
        #Font Sizes
        plt.setp(lines,linewidth=0.3)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.tight_layout()
    fig.savefig(plot_filename)
    sys.stdout.write("\r\t\t%i sequences processed (all finished)\n" % (count_seqs))
    sys.stdout.flush()

    #End program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("   PROGRAM FINISH TIME: " + time + '\n')

if __name__== "__main__":
    main()
