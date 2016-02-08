# -*- coding: utf-8 -*-
"""
Gene Finder Code!

@author: Liv Kelley

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

    I believe this is sufficient for DNA base pairs, but it would not cover RNA base pairs.
    If RNA is important it might be worth adding an if statement that paired Adenine with Uracil. 
    """
    # TODO: implement this

    if nucleotide  == 'A':
        return 'T'
    elif nucleotide  == 'C':
        return 'G'
    elif nucleotide  == 'G':
        return 'C'
    elif nucleotide  == 'T':
        return 'A'
    else:
        print 'Not a nitrogenous base.'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("")
    ''

    I think once again that these tests are sufficient if one is trying to find DNA base pairs.  
    The function of uracil is vital in producing proteins, though, and may be worth taking into account in some cases.
    """
    # TODO: implement this
    final_script = ''
    index = len(dna)-1
    while index > -1:
        letter = get_complement(dna[index])
        index = index - 1
        final_script = final_script+letter
    return final_script

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
    >>> rest_of_ORF('TGA') 
    ''
    >>> rest_of_ORF('')
    ''
    """

    for i in range(0, ((len(dna))/3)):
        n = 3*i 
        string = dna[n: n+3] 
        if string == 'TGA' or string == 'TAG' or string == "TAA": 
            return dna[0: n] 
    return dna  #This returns the whole string, because if there's no stop codons it shouldn't have to end it.  

def find_ATG (dna, search_point):
        n = search_point
        z = n+3 
        while z <= len(dna):
            if dna[n:z] == 'ATG':
                return n 
            else:
                n = n+3
                z = z+3 
        return -1

#Lauren worked with me on this problem, so if the first part of our code looks very similar, that is why.

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

    >>> find_all_ORFs_oneframe("CATCATATGCCTATGTAA")
    ['ATGCCTATG']

    >>> find_all_ORFs_oneframe('TGA')
    []

    >>> find_all_ORFs_oneframe('')
    []
    """  
    ORFs = []
    search_point = 0 
    n = find_ATG (dna, search_point)
    while n > -1: 
        ORFs.append(rest_of_ORF(dna[n:]))
        search_point = search_point + len(ORFs[-1]) + n 
        n = find_ATG (dna, search_point)
    return ORFs 

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

    >>> find_all_ORFs('TGATCATAG')
    []

    >>> find_all_ORFs('AT')
    []
    """
    alist = []
    frame1 = find_all_ORFs_oneframe(dna)
    alist.extend(frame1)
    frame2 = find_all_ORFs_oneframe(dna[1:])
    alist.extend(frame2)
    frame3 = find_all_ORFs_oneframe(dna[2:])
    alist.extend(frame3) 
    return alist 

#print(find_all_ORFs('ATGCATGAATGTAG')) 

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']

    >>> find_all_ORFs_both_strands('')
    []
    """
    blist = [] 
    blist.extend(find_all_ORFs(dna)) 
    blist.extend(find_all_ORFs(get_reverse_complement(dna)))
    return blist

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("")
    'No'
    >>> longest_ORF("AT")
    'No'
    """

    list_a = (find_all_ORFs_both_strands(dna))
    long_ORF = 0
    if list_a == []: 
        return "No"

    for i in range(0, len(list_a)):
        if len(list_a[i]) > len(list_a[long_ORF]):
            long_ORF = i

    return list_a[long_ORF]

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF

        >>> longest_ORF_noncoding('', 0)
        0

         """
    # TODO: implement this
    longer_ORF = 0 
    for i in range(0, num_trials):
        b = shuffle_string(dna)
        c = len(longest_ORF(b))
        if c >= longer_ORF: 
            longer_ORF = c
    return longer_ORF

def coding_strand_to_AA(dna): 

    """
        codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA('ATG')
        M 
        >>> coding_strand_to_AA('TAG')
        |
    """
 
    amino_acid_string = '' #Creates a list for compiling amino acids. 
    
    for i in range(0, len(dna)/3): #For the given codon
        codon = dna[3*i:3*i+3]
        amino_acid = aa_table[codon] #Convert to amino acids.
        amino_acid_string = amino_acid_string + amino_acid #Then add the result to the amino acid sequence.

    return amino_acid_string #The final amino acid sequence is the outcome.

#print(coding_strand_to_AA("TAG"))
       
def gene_finder(dna):

    """ Returns the amino acid sequences that are likely coded by the specified dna
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        >>> gene_finder('ATG')
        'M'

        >>> gene_finder('TAG')
        []

    """
    threshold = longest_ORF_noncoding (dna, 1500)
    coding_dna = list(find_all_ORFs_both_strands(dna)) 
    protein = []
    for i in coding_dna:
        if len(i) >= threshold: 
            protein.append(coding_strand_to_AA(i))
    return protein

dna = load_seq("./data/X73525.fa")
#print (gene_finder('TAG'))

if __name__ == "__main__":
     import doctest
     doctest.testmod() 



