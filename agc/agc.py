#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw
import operator

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    if (isfile(amplicon_file)):
        seq = ''
        for row in gzip.open(amplicon_file, 'r'):
            if row.startswith(">"):
                if len(seq)>= minseqlen:
                    yield seq
                seq = ''
            else:
                seq = seq + row.replace(" ", "").replace("\n", "") 
        yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    dico = {}
    for seq in read_fasta(amplicon_file, minseqlen):
            if seq in dico :
                dico[seq] += 1
            else:
                dico.update({seq: 1})
    
    sorted_d = sorted(dico.items(), key=operator.itemgetter(1),reverse=True)
    for i,j in (sorted_d):
        if j >= mincount:
            yield[i,j]

def get_chunks(sequence, chunk_size):
    chunks_list = []
    for i in range(0, len(sequence), chunk_size):
        if len(sequence[i:i+chunk_size]) == chunk_size:
            chunks_list.append(sequence[i:i+chunk_size])
        

    try :
        len(chunks_list) >= 4
    except : 
        raise ValueError
    else:
        return chunks_list 
    

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    for i in range(len(sequence) - kmer_size+1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    count = 0
    a, b = alignment_list
    for i in range(len(a)):
        if a[i] == b[i]:
            count += 1
    id = float(count)/float(len(a)) * 100
    return id  

def search_mates(kmer_dict, sequence, kmer_size):
    l = []
    res = [] 
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict.keys():
            for i in kmer_dict[kmer]:
                l.append(i)

    cnt = Counter(l).most_common(8)
    for i,j in cnt:
        res.append(i)
    return res

def detect_chimera(perc_identity_matrix):
    e = 0
    se1 = []
    se2 = []

    for p in perc_identity_matrix:
        e += statistics.stdev(p)
        se1.append(p[0])
        se2.append(p[1])

    if len(set(se1)) >= 2 or len(set(se2)) >= 2:
        ep = e/len(perc_identity_matrix)
        if ep > 5:
            return True
    return False

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    for i in cut_kmer(sequence, kmer_size):
        if i not in kmer_dict:
            kmer_dict[i] = list()
        kmer_dict[i].append(id_seq)
    return kmer_dict


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    file = open(output_file, "w")
    for i in range(len(OTU_list)):
        fasta = fill(OTU_list[i][0])
        file.write(">OTU_%d occurrence:%d\n" %(i+1, OTU_list[i][1]))
        file.write(fasta)
        file.write("\n")
    file.close()

#==============================================================
# Main program
#==============================================================

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    file = args.amplicon_file
if __name__ == '__main__':
    main()