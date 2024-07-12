from Bio import SeqIO

#initiializes lists of each neuc then creates the profile. byadding one in each spot where that neuc is detected
def profiler(seq_dict):
    first_key = next(iter(seq_dict))
    first_value = seq_dict[first_key]

    a = [0] * len(first_value)
    c = [0] * len(first_value)
    g = [0] * len(first_value)
    t = [0] * len(first_value)
    for key in seq_dict:
        for neuc in range(len(seq_dict[key])):
            if seq_dict[key][neuc] == 'A':
                a[neuc] += 1
            elif seq_dict[key][neuc] == 'C':
                c[neuc] += 1
            elif seq_dict[key][neuc] == 'G':
                g[neuc] += 1
            elif seq_dict[key][neuc] == 'T':
                t[neuc] += 1
    
    return a, c, g, t

#creates the consensus sequence from the max value of each neuc list for that partictular index. Has preference for ties being A>C>G>T
def concensusinator(a, c, g, t):
    consensus = []
    for i in range(len(a)):
        if a[i] == max(a[i], c[i], g[i], t[i]):
            consensus.append("A")
        elif c[i] == max(a[i], c[i], g[i], t[i]):
            consensus.append("C")
        elif g[i] == max(a[i], c[i], g[i], t[i]):
            consensus.append("G")
        elif t[i] == max(a[i], c[i], g[i], t[i]):
            consensus.append("T")

    output = ''.join(consensus)
    
    return output


if __name__ == "__main__":
   
   #this is biopython for 'Please take these FASTAs and make them into simple dicts, with the sequenece keyed by its name after the ">"'
    with open('data/cons.fas', 'r') as fasta_file:
        seq_records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
        seq_dict = {id: str(record.seq) for id, record in seq_records.items()}

    a, c, g, t = profiler(seq_dict)
    consensus = concensusinator(a, c, g, t)
    
    #Everything from here is formatting of the outputs
    a_strs = map(str, a)
    a_out = ' '.join(a_strs)
    c_strs = map(str,c)
    c_out = ' '.join(c_strs)
    g_strs = map(str, g)
    g_out = ' '.join(g_strs)
    t_strs = map(str, t)
    t_out = ' '.join(t_strs)


    print(consensus)
    print(f'A: {a_out}')
    print(f'C: {c_out}')
    print(f'G: {g_out}')
    print(f'T: {t_out}')






