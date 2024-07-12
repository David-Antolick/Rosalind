
from scipy import stats as forfunzies

#uses a binomial drstribution to find probability of at least n organisms in the 2^k generation being AaBb genotype
def bino_distribution(k, n):

    pop = 2 ** k
    prob_AaBb = 0.25
    out_prob = 1 - forfunzies.binom.cdf(n-1, pop, prob_AaBb)
    return out_prob



if __name__ == "__main__":
    with open('data/lia.txt', 'r') as raw:
        numbs = raw.readline().split(' ')
        k, n = map(int, numbs)

    out = bino_distribution(k, n)
    
    #formatting and output
    print (f'This is K:{k} and this is N ():{n}')
    print(round(out, 3))