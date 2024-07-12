from Bio.SeqUtils import molecular_weight

#Calculates total protein weight from monoisometic weights (from ROSALIND) and returns rounded to thousandths
def calc_monoiso_mass(prtn_str):
    masses = {
        'A':71.03711,
        'C':103.00919,
        'D':115.02694,
        'E':129.04259,
        'F':147.06841,
        'G':57.02146,
        'H':137.05891,
        'I':113.08406,
        'K':128.09496,
        'L':113.08406,
        'M':131.04049,
        'N':114.04293,
        'P':97.05276,
        'Q':128.05858,
        'R':156.10111,
        'S':87.03203,
        'T':101.04768,
        'V':99.06841,
        'W':186.07931,
        'Y':163.06333
    }

    weight = 0
    for let in prtn_str:
        weight += masses[let]
    round_weight = round(weight, 3)
    return round_weight


if __name__ == '__main__':
    with open('data/prtm.txt', 'r') as raw:
        prtn_str = raw.readline().strip().upper()
 
    print(calc_monoiso_mass(prtn_str))

    #this uses the monosiotopic weights in the bio library, yet gives the wrong answer. Probably difference in measurements
    '''prtn_weight = molecular_weight(prtn_str, 'protein', monoisotopic=True)
    print(round(prtn_weight, 2))'''


    