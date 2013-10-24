#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def usage():
    test="name"
    message='''
python ChunkGenome.py --reference ../input/MSU_r7.fa --assembly assembly.fa --project HEG4.ALLPATHLG
Cut the contig of assembly into 1kb pieces and mapped onto reference genome to evaluate the misassembly rate of assembly.

    '''
    print message

'''convert assembly into gap-free contigs, contigs less than 1 kb will not output and not used for evalution'''
def assembly2contig(assembly, output):
    if os.path.exists(output):
        print output + ' exists, delete it or rename the output dirctory'
        sys.exit(2)
    else:
        os.system('mkdir ' + output)

    with open (output + '/assembly.contig.fa', 'w') as ofile:
        for record in SeqIO.parse(assembly, 'fasta'):
            contigs = re.split(r'N+',str(record.seq))

            for i in range(len(contigs)):
                ctgid = record.id + '_' + str(i)
                ctgseq= Seq(contigs[i])
                ctg   = SeqRecord(ctgseq, id=ctgid, description='')
                if len(contigs[i]) >= 1000:
                    SeqIO.write(ctg, ofile, 'fasta')

'''chunk the contig(Seq object) into 1kb sequence (Seq object), chunks smaller than 1kb will be discard'''
def chunk(contig, ctgid):
    #contigstr = str(contig)
    chunks = defaultdict(list)
    win = 1000
    run = int(len(str(contig))/win)
    for i in range(run):
        start = i * win
        end   = (i + 1) * win - 1
        chk   = contig[start:end]
        #print ctgid, run, start, end, '\n', chk
        chkid = ctgid + '_chk' + str(i)
        chunks[chkid] = chk
    return chunks

def chunk_genome(output):
    contig = output + '/assembly.contig.fa'
    with open (output + '/assembly.contig.1k_chunk.fa', 'w') as ofile:
        for record in SeqIO.parse(contig, 'fasta'):
            #print record.id, record.seq
            chunks = chunk(record.seq, record.id)
            for ctgid in sorted(chunks.keys()):
                chk = SeqRecord(chunks[ctgid],id=ctgid, description='')
                SeqIO.write(chk, ofile, 'fasta')

def map_sequence(reference, output):
     


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference')
    parser.add_argument('-a', '--assembly')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.assembly) > 0 and len(args.reference) > 0
    except:
        usage()
        sys.exit(2)

    if args.project is None:
        args.project = 'Evaluation_Chunk'

    assembly2contig(args.assembly, args.project)
    chunk_genome(args.project)
    map_sequence(args.reference, args.project)
    #summary(args.project)

if __name__ == '__main__':
    main()

