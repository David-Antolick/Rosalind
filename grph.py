from Bio import SeqIO

#Addes 2 more values to each key being the prefix and suffix of length 3
def _dict_organizer(fasta_dict):
    new_dict = {}

    k = 3
    for key in fasta_dict:
        prefx = fasta_dict[key][:k]
        suffx = fasta_dict[key][-k:]
        new_dict[key] = [fasta_dict[key], prefx, suffx]
    
    return new_dict

def adjacency_detector(fasta_dict):
    
    dict_plus = _dict_organizer(fasta_dict)
    adj_list = []

    for key1 in dict_plus:
        for key2 in dict_plus:
            if key2 != key1:
                if dict_plus[key1][1] == dict_plus[key2][2]:
                    adj_list.append([key1, key2])

    return adj_list




if __name__ == "__main__":
   
   #this is biopython for 'Please take these FASTAs and make them into simple dicts, with the sequenece keyed by its name after the ">"'
    with open('data/grph.fas', 'r') as fasta_file:
        seq_records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
        seq_dict = {id: str(record.seq) for id, record in seq_records.items()}
    print(seq_dict)

    adj_list = adjacency_detector(seq_dict)
    #printing made neat (the suffix containing sequence must be printed first or the answer will not be accepted)
    for list in adj_list:
        print(list[1], list[0])
