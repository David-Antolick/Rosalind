

def mrna(input):
    codons = {
                'F': ['UUU', 'UUC'],
                'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
                'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                'Y': ['UAU', 'UAC'],
                '*': ['UAA', 'UAG', 'UGA'],
                'C': ['UGU', 'UGC'],
                'W': ['UGG'],
                'P': ['CCU', 'CCC', 'CCA', 'CCG'],
                'H': ['CAU', 'CAC'],
                'Q': ['CAA', 'CAG'],
                'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                'V': ['GUU', 'GUC', 'GUA', 'GUG'],
                'A': ['GCU', 'GCC', 'GCA', 'GCG'],
                'D': ['GAU', 'GAC'],
                'E': ['GAA', 'GAG'],
                'G': ['GGU', 'GGC', 'GGA', 'GGG'],
                'I': ['AUU', 'AUC', 'AUA'],
                'M': ['AUG'],
                'T': ['ACU', 'ACC', 'ACA', 'ACG'],
                'N': ['AAU', 'AAC'],
                'K': ['AAA', 'AAG']
            }

    cnt = 1
    for aa in protein:
        cnt = cnt * len(codons[aa])
    cnt*=3
    out = cnt % 1000000
    return out










if __name__ == "__main__":
    f = open('data/mrna.txt', 'r')
    protein = f.readline().strip()
    out = mrna(protein)
    print(out)