import reuseable as RU
import re

#this returns a dict containing the proteinID: motif positions. It only returns IDs with matching motifs.
#this function will NOT return overlapping motifs
def N_glyco_motif(prtn_dict):
    motif = re.compile(r'N[^P][ST][^P]')
    matches_dict = {}
    for id in prtn_dict:
        strt_pos = []
        for aa in range(len(prtn_dict[id])):
            if motif.match(prtn_dict[id], aa):
                strt_pos.append(aa+1)
            if strt_pos:
                matches_dict[id] = strt_pos
    return matches_dict
    


if __name__ == "__main__":
    with open('data/mprt.fas', 'r') as raw:
        prtn_list = []
        for prtn in raw.readlines():
                cleasned = prtn.strip().split('_')[0]
                prtn_list.append(cleasned)

    #fuction calls
    prtn_dict = RU.uniprot_todict(prtn_list, 'product/mprt_write.fas')
    matches_dict = N_glyco_motif(prtn_dict)

    #Formatting of outputs
    #Doesn't include additional text on some answers, but it has whats important
    for id in matches_dict:
        str_list = [str(i) for i in matches_dict[id]]
        print(f'{id} \n{" ".join(str_list)}')