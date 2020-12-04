import statistics
import math
import numpy as np


class Population():
    """Model a biological population."""

    def __init__(self, pfreq, heterozygosity=None, AAfreq=None, Aafreq=None, twoN=100, mutation_rate=None,
                 wAA=None, wAa=None, waa=None):
        """Set allele and genotype frequencies."""
        self.pfreq = pfreq
        self.twoN = twoN
        if heterozygosity is not None:
            self.Aafreq = heterozygosity
        if AAfreq and Aafreq:
            self.AAfreq = AAfreq
            self.Aafreq = Aafreq
            self.aafreq = 1-AAfreq-Aafreq
        self.wAA = wAA
        self.wAa = wAa
        self.waa = waa
        if mutation_rate:
            self.theta = 4*(twoN/2)*mutation_rate

    @property
    def qfreq(self):
        return round(1-self.pfreq, 2)

    def get_hw_freqs(self):
        """Get Hardy-Weinberg frequencies."""
        hwFreqs = {'AA': round((self.pfreq)**2, 2),
                   'Aa': round(2*self.pfreq*self.qfreq, 2),
                   'aa': round((self.qfreq)**2, 2)}
        return hwFreqs

    def get_fixation_index(self):
        """Get fixation index."""
        expectedH = 2*self.pfreq*self.qfreq
        observedH = self.Aafreq
        fixation_index = (expectedH-observedH)/expectedH
        return round(fixation_index, 2)

    def theta_hetero(self):
        """Calculate expected heterozygosity based on theta."""
        thetaHetero = self.theta/(self.theta + 1)
        return thetaHetero

    def drift_sel_balance(self, selection):
        """See drift-selection balance."""
        fourNeS = 4*self.twoN*selection
        if fourNeS > 1:
            return f'4Nes: {fourNeS}\n' \
                'Selection is stronger than drift.'
        elif fourNeS < 1:
            return f'4Nes: {fourNeS}\n' \
                'Drift is stronger than selection.'
        else:
            return f'4Nes: {fourNeS}\n' \
                'Drift and selection are equal.'

    def prob_coalescence_discrete(self, t):
        """Obtain the probability of coalescence given t discrete time intervals."""
        prob_not = (1-(1/self.twoN))**(t-1)
        prob_do = 1/self.twoN
        prob = prob_not*prob_do
        return round(prob, 4)

    def sim_drift(self, gens):
        """Simulate genetic drift, using a binomial distribution."""
        init_freq = self.pfreq
        freqs_vector = [init_freq]
        for generation in np.arange(gens):
            if self.pfreq >= 0.995 or self.pfreq <= 0.005:
                freqs_vector.append(round(self.pfreq))
            else:
                new_pfreq = (np.random.binomial(self.twoN, self.pfreq))/self.twoN
                if new_pfreq > 1 or new_pfreq < 0:
                    freqs_vector.append(round(self.pfreq))
                    self.pfreq = new_pfreq
                else:
                    freqs_vector.append(new_pfreq)
                    self.pfreq = new_pfreq
        self.pfreq = round(freqs_vector[-1], 4)


    def get_mix_stats(self, otherPop):
        """Get F_IS, F_ST, and F_IT for two populations."""
        h_I1 = self.Aafreq
        h_I2 = otherPop.Aafreq
        h_I = (h_I1 + h_I2)/2

        h_S1 = self.get_hw_freqs()['Aa']
        h_S2 = otherPop.get_hw_freqs()['Aa']
        h_S = (h_S1 + h_S2)/2

        p_bar = (self.pfreq + otherPop.pfreq)/2
        q_bar = 1-p_bar
        h_T = 2*p_bar*q_bar

        structure = {'F_IS': round((h_S - h_I)/h_S, 2),
                    'F_ST': round((h_T - h_S)/h_T, 2),
                    'F_IT': round((h_T - h_I)/h_T, 2)}

        return structure

    def drift_geneFlow_balance(self, m):
        """See drift-geneFlow balance."""
        fourNem = 4*self.twoN*m
        if fourNem > 1:
            return f'4Nem: {fourNem}\n' \
                'Gene flow is stronger than drift.'
        elif fourNem < 1:
            return f'4Nem: {fourNem}\n' \
                'Drift is stronger than gene flow.'
        else:
            return f'4Nem: {fourNem}\n' \
                'Drift and gene flow are equal.'

    def future_hetero(self, t):
        """Get heterozygosity after one time period."""
        newHetero = ((1-(1/self.twoN))**t)*self.Aafreq
        return round(newHetero, 4)

    def sim_selection(self):
        """Simulate 1 generation of natural selection."""
        numerator = (((self.pfreq)**2)*self.wAA) + ((self.pfreq*self.qfreq)*self.wAa)
        denominator = (((self.pfreq)**2)*self.wAA) + (2*(self.pfreq*self.qfreq)*self.wAa) + (((self.qfreq)**2)*self.waa)
        pt2 = numerator/denominator
        delta_p = pt2-self.pfreq
        self.pfreq = round(pt2, 4)
        return f'Delta p: {delta_p:.4f}'


def get_diseq_stats(g11, g12, g21, g22, p1=None, p2=None):
    """Get gametic disequilibrium coefficient, D, as well as r^2."""
    D = (g11*g22)-(g12*g21)
    if p1 and p2:
        q1 = 1-p1
        q2 = 1-p2
        r_2 = D/np.sqrt(p1*q1*p2*q2)
        return f'D: {D:.4f}\n' \
            f'r2: {r_2:.4f}'
    else:
        return f'D: {D:.4f}'
   
def Ne_through_time(*args):
    """Obtain estimate of Ne by."""
    Ne = statistics.harmonic_mean(*args)
    return Ne

def solve_breeders(R=None, h2=None, S=None):
    """Solve the Breeder's equation given two variables."""
    if R and h2:
        S = R/h2
        return f'S: {S:.2f}'
    elif R and S:
        h2 = R/S
        return f'h2: {h2:.2f}'
    else:
        R = h2*S
        return f'R: {R:.2f}'

def get_mean_fitness(pfreq, wAA, wAa, waa):
    """Get mean fitness."""
    meanFitness = (((pfreq)**2)*wAA) + (2*(pfreq*(1-pfreq)*wAa)) + (((1-pfreq)**2)*waa)
    return round(meanFitness, 4)

def exp_growth(n0, r, t):
    """Get future size of population after exponential growth."""
    futureSize = n0*math.exp(r*t)
    return round(futureSize)

def inf_island_fst(Ne, m=None, t=None):
    """Get Fst in the infinite island model, given m or t."""
    if m:
        Fst = 1/((4*Ne*m)+1)
        return round(Fst, 4)
    else:
        Fst = 1-(math.exp(-t/2*Ne))
        return round(Fst, 4)

def inf_island_Nem(Fst):
    """Get Nem in the infinite island model, given Fst."""
    Nem = (1/4)*((1/Fst)-1)
    return round(Nem, 2)

def colony_fst(Fst, k, phi):
    """Get colony Fst.""" 
    colonyFst = 1/(2*k)+(phi*(1-(1/(2*k)))*Fst)
    return round(colonyFst, 4)


def genetic_assignment(pop1, pop2, pop3, genotype):
    """
    Each pop should be a dictionary with keys 'A', 'a', 'B', 'b', and 'Prior.'
    They should *not* be population objects from the Population class in this module.
    
    Genotype should be specified as a list of alleles, e.g., ['A', 'A', 'B', 'B'].
    """
    locus1genotype  = genotype[0:2]
    locus2genotype = genotype[2:4]

    allelesFreqsPop1 = {
        'A': pop1['A'],
        'a': pop1['a'],
        'B': pop1['B'],
        'b': pop1['b']
    }

    locus1conditionalPop1 = allelesFreqsPop1[locus1genotype[0]]*allelesFreqsPop1[locus1genotype[1]]
    locus2conditionalPop1 = allelesFreqsPop1[locus2genotype[0]]*allelesFreqsPop1[locus2genotype[1]]
    multilocusConditionalPop1 = locus1conditionalPop1*locus2conditionalPop1

    locus1genotypeProbPop1 = locus1conditionalPop1*pop1['Prior']
    locus2genotypeProbPop1 = locus2conditionalPop1*pop1['Prior']
    multilocusGenotypeProbPop1 = multilocusConditionalPop1*pop1['Prior']

    allelesFreqsPop2 = {
        'A': pop2['A'],
        'a': pop2['a'],
        'B': pop2['B'],
        'b': pop2['b']
    }

    locus1conditionalPop2 = allelesFreqsPop2[locus1genotype[0]]*allelesFreqsPop2[locus1genotype[1]]
    locus2conditionalPop2 = allelesFreqsPop2[locus2genotype[0]]*allelesFreqsPop2[locus2genotype[1]]
    multilocusConditionalPop2 = locus1conditionalPop2*locus2conditionalPop2

    locus1genotypeProbPop2 = locus1conditionalPop2*pop2['Prior']
    locus2genotypeProbPop2 = locus2conditionalPop2*pop2['Prior']
    multilocusGenotypeProbPop2 = multilocusConditionalPop2*pop2['Prior']

    allelesFreqsPop3 = {
        'A': pop3['A'],
        'a': pop3['a'],
        'B': pop3['B'],
        'b': pop3['b']
    }

    locus1conditionalPop3 = allelesFreqsPop3[locus1genotype[0]]*allelesFreqsPop3[locus1genotype[1]]
    locus2conditionalPop3 = allelesFreqsPop3[locus2genotype[0]]*allelesFreqsPop3[locus2genotype[1]]
    multilocusConditionalPop3 = locus1conditionalPop3*locus2conditionalPop3

    locus1genotypeProbPop3 = locus1conditionalPop3*pop3['Prior']
    locus2genotypeProbPop3 = locus2conditionalPop3*pop3['Prior']
    multilocusGenotypeProbPop3 = multilocusConditionalPop3*pop3['Prior']

    locus1posteriorProbs = {
        'Pop1': round((locus1conditionalPop1*pop1['Prior']) / 
        (locus1genotypeProbPop1+locus1genotypeProbPop2+locus1genotypeProbPop3), 2),
        'Pop2': round((locus1conditionalPop2*pop2['Prior']) / 
        (locus1genotypeProbPop1+locus1genotypeProbPop2+locus1genotypeProbPop3), 2),
        'Pop3': round((locus1conditionalPop3*pop3['Prior']) / 
        (locus1genotypeProbPop1+locus1genotypeProbPop2+locus1genotypeProbPop3), 2)
    }

    locus2posteriorProbs = {
        'Pop1': round((locus2conditionalPop1*pop1['Prior']) / 
        (locus2genotypeProbPop1+locus2genotypeProbPop2+locus2genotypeProbPop3), 2),
        'Pop2': round((locus2conditionalPop2*pop2['Prior']) / 
        (locus2genotypeProbPop1+locus2genotypeProbPop2+locus2genotypeProbPop3), 2),
        'Pop3': round((locus2conditionalPop3*pop3['Prior']) / 
        (locus2genotypeProbPop1+locus2genotypeProbPop2+locus2genotypeProbPop3), 2)
    }

    multilocusPosteriorProbs = {
        'Pop1': round((multilocusConditionalPop1*pop1['Prior']) / 
        (multilocusGenotypeProbPop1+multilocusGenotypeProbPop2+multilocusGenotypeProbPop3), 2),
        'Pop2': round((multilocusConditionalPop2*pop2['Prior']) / 
        (multilocusGenotypeProbPop1+multilocusGenotypeProbPop2+multilocusGenotypeProbPop3), 2),
        'Pop3': round((multilocusConditionalPop3*pop3['Prior']) / 
        (multilocusGenotypeProbPop1+multilocusGenotypeProbPop2+multilocusGenotypeProbPop3), 2)
    }

    return f'A Locus: {locus1posteriorProbs}\nB Locus: {locus2posteriorProbs}\nMultilocus: {multilocusPosteriorProbs}'
    

