# -*- coding: utf-8 -*-
"""
Gene Finder mini project for Olin Software Design Spring 2016
The base of this code was not written by the author

Date last worked on: 1/19/16

@author: Ian Paul

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    #TODO use the list index trick
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO:
    complement = ''
    for i in dna:
        complement = get_complement(i) + complement
    return complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: use less data, be faster. find a functional way to do this.
    stop_codons = ['TAA', 'TAG', 'TGA'] #should probably be defined globally
    for i in range(0,len(dna),3):
        codon = dna[i:i+3]
        for stop in stop_codons:
            if (codon == stop):
                return dna[:i]
    return dna



def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this
    stop_codons = ['TAA', 'TAG', 'TGA'] #should probably be defined globally
    start_codon = 'ATG'
    orfs = []
    i =0

    while i < len(dna):
        codon = dna[i:i+3]
        if (codon ==start_codon):
            orf = rest_of_ORF(dna[i:])
            orfs.append(orf)
            i = i + len(orf)
        i = i +3
    return orfs





def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this
    return [i for j in range(0,3) for i in find_all_ORFs_oneframe(dna[j:])]
    #list comprehension is cool


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    longest = ''
    for orf in find_all_ORFs_both_strands(dna):
    	if len(orf) > len(longest):
    		longest = orf
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
         """
    # TODO: find some other way to test, tested by hand, seems legit
    res = 0
    for i in range(num_trials):
    	working_dna = shuffle_string(dna)
    	if len(longest_ORF(working_dna)) > res:
    		res = len(longest_ORF(working_dna))
    return res



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    acid = ''
    for i in range(0,len(dna)/3):
    	acid = acid + aa_table[dna[3*i:3*(i+1)]]
    return acid


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    threshold = longest_ORF_noncoding(dna,1500)
    orfs = find_all_ORFs_both_strands(dna)
    genes = []
    for orf in orfs:
    	if orf < threshold:
    		orfs.remove(orf)
    for orf in orfs:
    	genes.append(coding_strand_to_AA(orf))
    return genes


if __name__ == "__main__":
    import doctest
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print(dna)
    doctest.testmod()
    #print(gene_finder(dna))