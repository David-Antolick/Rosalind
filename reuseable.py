import requests
from Bio import SeqIO
from Bio.Seq import Seq
import re


#the next two functions are able to create a simplified dict from a uniprot text
#They first write to a text file then convert to an id:sequence style dict
def _extractor(id):
     return id.split('|')[1]

def uniprot_todict(prtn_list, write_file):
    #'http://www.uniprot.org/uniprot/*PROTEIN*.fasta.'

    with open(write_file, "w") as f:
        for prtn in prtn_list:
            response = requests.get(f'https://rest.uniprot.org/uniprotkb/{prtn}.fasta')
            assert response.status_code == 200, f'API returned {response.status_code} \n {response.text}'
            info = response.text
            f.write(f'{info}\n')

    #Creates a dict in the stle of prtnID:sequence from write_file raw data
    with open(write_file, 'r') as fasta_file:
        seq_records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
        seq_dict = {_extractor(id): str(record.seq) for id, record in seq_records.items()}

    print('uniprot_todict: File created!')
    return seq_dict




#converts RNA to DNA one string at a time
def RNAtoDNA(seq):
    seq_uppies = seq.upper().strip()
    seq_rna = re.sub('U', 'T', seq_uppies)
    return seq_rna

#converts DNA to RNA one string at a time
def DNAtoRNA(seq):
    seq_uppies = seq.upper().strip()
    seq_rna = re.sub('T', 'U', seq_uppies)
    return seq_rna
    
    
#reads RNA to identify ORFs, returns a dict of ORFs by ID, As well as a dict by protein

'''def readingframes(rna_dict):

    pattern = re.compile(r'AUG(?:[ACGU]{3})*?(UAG|UAA|UGA)')

    orf_dict_rna = {}
    for key in rna_dict:
        orfs = [match.group() for match in pattern.finditer(rna_dict[key])]
        orf_dict_rna[key] = orfs
    
    orf_dict_prtn = {}
    for key in orf_dict_rna:
        prtn_lst = []
        for orf in orf_dict_rna[key]:
            rna = Seq(orf)
            translation = str(rna.translate())
            prtn_lst.append(translation)
        orf_dict_prtn[key] = prtn_lst
    
    

    
    return orf_dict_prtn'''

#takes an rna dictionary, returns list of DISTINCT proteins from it (and its reverse comliment, ofc)
def openreadingframes(rna_dict):
    

    #creates a protein dictionary with the three possible translations (biopython is upset about this, I think it cries every time its run)
    prtn_dict = {}
    translations = []

    #translates to AA codes in each reading frame
    for key in rna_dict:
        for i in range(3):
            rna = Seq(rna_dict[key][i:])
            translations.append(str(rna.translate()))
        prtn_dict[key] = translations
        translations = []

    #discovers ALL possible protein outputs following the regex patter r'M[^*]*\*',
    pattern = re.compile(r'M[^*]*\*')

    orf_dict = {}
    orf_list = []

    for key in prtn_dict:
        for i in range(3):
            orf_list.extend([match.group() for match in pattern.finditer(prtn_dict[key][i])])
        orf_dict[key] = orf_list
        print(key, prtn_dict[key])
        orf_list = []
    
    #checks to see if there are any proteins hiding within other proteins
    for key in orf_dict:
        for string in range(len(orf_dict[key])):
            for aa in range(len(orf_dict[key][string])):
                if aa != 0:
                    if orf_dict[key][string][aa] == 'M':
                            orf_dict[key].append(orf_dict[key][string][aa:])


    #merges the dictionary keys of the original and reverse compliment of the RNA together
    #still allows for multiple initial RNA strings to be handled
    #also picks out duplicates within the original and reverse compliment
    key_list = list(orf_dict.keys())
    merged_dict = {}

    for key in range(0, len(orf_dict), 2):
        merged_dict[key_list[key]] = list(set(orf_dict[key_list[key]] + orf_dict[key_list[key + 1]]))

    return merged_dict


    