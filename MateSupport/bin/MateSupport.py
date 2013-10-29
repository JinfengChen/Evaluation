#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
import pysam

def usage():
    test="name"
    message='''
python MateSupport.py --bam ../input/out.smaltmap.bam --type inward --insert 3000 --sd 400 --pair 12188423 
--bam: input bam file
--type: pair of read is inward or outward. read retrieved from allpathlg are all inward. DEFAULT is inward and we only deal with inward.
--insert: insertion size for library, use allpathlg estimated size
--sd: insertion size stardard deviation for library, use allpathlg estimated size 
--pair: total pair of read in fasta used to map, if pair not given we use paired in mapped reads.
--project: prefix of output
Use "cat /rhome/cjinfeng/BigData/00.RD/Evaluation/REAPR/input/HEG4/jump_reads_ec.illuminaGAIIx_3k.A.fastq | awk '{if(NR%4==2) print length($1)}'| wc -l" to get pair from fasta

Give a bam file of mating pair library, we calculate the number of pairs that mapped properly on the assembly.
    '''
    print message


def parsebam(bam, type, insert, sd, pair, prefix):
    support = 0
    wrong   = 0
    allpair = 0
    unpair  = 0
    total   = 0
    pool = defaultdict(int)
    minsize = 0 if insert - sd * 6 < 0 else insert - sd * 6
    maxsize = insert + sd * 6
    samfile = pysam.Samfile(bam, 'rb')
    '''fetch only get mapped read in bam file'''
    for read in samfile.fetch():
        total += 1
        '''skip at second read in pair and delete this key from hash table as we only have one pair of this reads'''
        if pool.has_key(read.qname):
            del pool[read.qname]
            continue
        
        if read.is_paired:
            read1_name   = read.qname
            read1_refid  = read.tid
            read1_pos    = read.pos
            read1_strand = 1 if read.is_reverse else 0
            read2_refid  = read.rnext
            read2_pos    = read.pnext
            read2_strand = 1 if read.mate_is_reverse else 0
            '''store read name in hash table when meet the first reads in pair'''
            pool[read1_name] =1
            '''as the sam file are sorted by coordinate, we always expected the first come on the forward strand'''
            if read1_refid == read2_refid and read1_strand == 0 and read2_strand == 1 and read.tlen >= minsize and read.tlen <= maxsize:
                support += 1
                #print '>properly pair'
                #print read
                #print read1_name, read1_refid, read1_pos, read1_strand
                #print read1_name, read2_refid, read2_pos, read2_strand
                #print read.tlen
            else:
                wrong += 1
                #print '>not properly pair'
                #print read
                #print read1_name, read1_refid, read1_pos, read1_strand
                #print read1_name, read2_refid, read2_pos, read2_strand
                #print read.tlen
        else:
            unpair += 1   

    samfile.close()
    totalpair = support + wrong if pair is None else pair
    ofile = open (prefix + '.sum', 'w')
    if pair:
        print >> ofile, 'Total pairs in fasta:', totalpair
        print >> ofile, 'Supporting pairs: ', support, '%.2f' % (100.00 * float(support)/float(pair))
        print >> ofile, 'Not properly mapped: ', wrong, '%.2f' % (100.00 * float(wrong)/float(pair))
        print >> ofile, 'Total pairs mapped: ', support + wrong
        print >> ofile, 'Supporting pairs: ', support, '%.2f' % (100.00 * float(support)/float(support + wrong)) 
        print >> ofile, 'Not properly mapped: ', wrong, '%.2f' % (100.00 * float(wrong)/float(support + wrong))
    else:
        print >> ofile, 'Total pairs mapped: ', support + wrong
        print >> ofile, 'Supporting pairs: ', support, '%.2f' % (100.00 * float(support)/float(support + wrong))
        print >> ofile, 'Not properly mapped: ', wrong, '%.2f' % (100.00 * float(wrong)/float(support + wrong))

    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam')
    parser.add_argument('-t', '--type')
    parser.add_argument('-i', '--insert', type=int)
    parser.add_argument('-s', '--sd', type=int)
    parser.add_argument('-p', '--pair', type=int)
    parser.add_argument('-P', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bam) > 0
    except:
        usage()
        sys.exit(2)

    if args.type is None:
        args.type = 'inward'
    if args.insert is None or args.sd is None:
        usage()
        sys.exit()
    

    parsebam(args.bam, args.type, args.insert, args.sd, args.pair, args.project)

if __name__ == '__main__':
    main()

