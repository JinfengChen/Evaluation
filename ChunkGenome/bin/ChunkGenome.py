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

Example output:
class  errorrate count  percent(%)
1	0%	59432	76.98
2	<= 1%	13876	17.97
3	<= 10%	2901	3.76
4	> 10%	760	0.98
5	No match	235	0.30

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
    cmdline = 'perl runblat.pl --infile ' + output + '/assembly.contig.1k_chunk.fa --database ' + reference 
    os.system(cmdline)
    os.system('mv assembly.contig.1k_chunk.fa.best.psl ' + output)  
 
def fasta_num(fasta):
    num = 0
    for record in SeqIO.parse(fasta,"fasta"):
        num += 1
    return num
 

'''998     1       0       0       0       0       1       1       +       scaffold_0_0_chk11      999     0       999     Chr10   23207287        5199409 5200409 2       867,132,        0,867,  5199409,5200277,'''
def summary(output):
    sumdata = defaultdict(int)
    errorgp = {'1' : '0.0', '2' : '<= 1%', '3' : '<= 10%', '4' : '> 10%', '5' : 'No match'}
    allhit  = 0
    best = output + '/assembly.contig.1k_chunk.fa.best.psl'
    ofile = open(output + '/assembly.contig.1k_chunk.match.info', 'w')
    print >> ofile, 'Chunk\tLength\tMatch\tMismatch\tErrorRate\tClass\n' 
    with open (best, 'r') as infile:
       for line in infile:
           line = line.rstrip()
           unit = re.split(r'\t', line)
           allhit += 1
           error= float((float(unit[10])-float(unit[0]))/float(unit[10]))
           group= 0
           ''' == 0.0 is perfect; <=0.001 and <= 0.01 is match; <= 0.1 and <= 0.5 is misassembly; > 0.5 is novel sequence '''
           '''1 is perfect; 2 is match; 3 is misassembly; 4 is novel sequence'''
           if error == 0.0:
               group = 1
           elif error <= 0.001:
               group = 2
           elif error <= 0.01:
               group = 2
           elif error <= 0.1:
               group = 3
           elif error <= 0.5:
               group = 4
           else:
               group = 4
           sumdata[group] = sumdata[group] + 1
           out  = [unit[9], unit[10], unit[0], unit[1], error, group]
           out  = map(str, out)
           print >> ofile, '\t'.join(out)
    ofile.close()
    
    ofile = open(output + '/assembly.contig.1k_chunk.match.sum', 'w')
    allseq = fasta_num(output + '/assembly.contig.1k_chunk.fa')
    sumdata['5'] = allseq-allhit
    for g in sorted(sumdata.keys()):
        rate = 100.00 * float(sumdata[g])/float(allseq)
        out = [g, errorgp[str(g)], sumdata[g], "%.2f" % rate]
        out = map(str, out)
        print >> ofile, '\t'.join(out)
    ofile.close()
    #print allseq, allhit, allseq-allhit, (float(allseq)-float(allhit))/float(allseq)

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
    summary(args.project)

if __name__ == '__main__':
    main()

