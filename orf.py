from Bio import SeqIO
from Bio.Seq import Seq
import reuseable as RU




if __name__ == "__main__":
    with open('data/orf.fas', 'r') as fas_item:
        seq_records = SeqIO.to_dict(SeqIO.parse(fas_item, 'fasta'))
        seq_dict = {id: str(record.seq) for id, record in seq_records.items()}
    
    
    #convert seqs to RNA AND addes its reverse compliment
    rna_dict = {}
    for key in seq_dict:
        rna_dict[key] = RU.DNAtoRNA(seq_dict[key])
        rna_dict[f'{key}_r'] = RU.DNAtoRNA(str((Seq(seq_dict[key]).reverse_complement())))
    
    print(rna_dict)
    new_dict = RU.openreadingframes(rna_dict)

    for key in new_dict:
        for i in range(len(new_dict[key])):
            print(new_dict[key][i][:-1])

