from reuseable import basic as base
from reuseable import DNA as dna
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq

def rev_palindrome(seq_dict):
    search_range = range(4, 13)
    palin_list = []

    #primary function 
    for id in seq_dict:
        curr = seq_dict[id]

        for i in search_range:
            
            for aa in range(len(curr) - i + 1):
                seq = Seq(curr[aa:aa+i])
                rev_comp = seq.reverse_complement()

                if seq == rev_comp:
                    palin_list.append((aa+1, len(seq)))
    
    return palin_list


        

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
    seq_dict = base.fasta_simp_dict('data/revp.fas')


    palin_list = rev_palindrome(seq_dict)

    for toup in palin_list:
        print(f'{toup[0]} {toup[1]}')
