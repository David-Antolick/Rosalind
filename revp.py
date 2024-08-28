from reuseable import basic as base
from reuseable import DNA as dna
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq

def rev_palindrome(seq_dict):
    search_range = range(4, 13)
    palin_list = []

    #primary function 
    for id in seq_dict:
        fwd = seq_dict[id][0]
        bkw = seq_dict[id][1]
        leng = len(fwd)
        leng_bkw = len(bkw)
        if leng != leng_bkw:
            raise Exception('Length of sequence and reverse compliment do not match.')
        for i in search_range:
            for aa in range(len(fwd)):
                seq = Seq(fwd[aa:aa+i])
                if len(seq) == i:
                    result = nt_search(fwd, seq)
                    if result[1] != aa:
                        palin_list.append((len(result[0]), result[1]))
    print(palin_list)


        

    #this is a grade-A mess, no idea tbh
    '''for id in seq_dict:
        fwd = seq_dict[id][0]
        bkw = seq_dict[id][1]
        leng = len(fwd)
        leng_bkw = len(bkw)
        if leng != leng_bkw:
            raise Exception('Length of sequence and reverse compliment do not match.')
        for i in search_range:
            for aa in range(len(fwd)-i):
                if fwd[aa] == bkw[aa + i]:
                    if fwd[aa+1] == bkw[aa + i - 1]:
                        print(aa)
                else: pass'''
                    
                    




if __name__ == "__main__":
    fwd_dict = base.fasta_simp_dict('data/revp.fas')

    seq_dict = dna.reverse_complement(fwd_dict)

    rev_palindrome(seq_dict)
