from Bio import SeqIO
import json
import pandas as pd

with open('cons.fas', 'r') as fasta_file:
    seq_records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    seq_dict = {id: str(record.seq) for id, record in seq_records.items()}

print(json.dumps(seq_dict, indent = 2))
"""transfrom into one column, use dummys, create n spreadsheets for each thing, vert stack spreadsheets, then find what column is the maX VALUE FOR EACH ROW TO CREATE CONSENSUS"""