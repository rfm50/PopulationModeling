import numpy as np
import matplotlib.pyplot as plt


def drift_simulator(NeSmall=10, NeMedium=100, NeLarge=1000, replicates=10,
                    gens=100, startingFreq=0.5):

    fig, axs = plt.subplots(1, 3, figsize=(20, 5))

    for ind, twoN in enumerate([2*NeSmall, 2*NeMedium, 2*NeLarge]):
        

        pops = []
        for replicate in np.arange(replicates):
            p = startingFreq
            q = 1-p

            freqs = [p]

            for generation in np.arange(gens):
                if p >= 0.995 or p <= 0.005:
                    freqs.append(round(p))
                else:
                    newP = (np.random.binomial(twoN, p))/twoN
                    if newP > 1 or newP < 0:
                        freqs.append(round(p))
                        p = newP
                    else:
                        freqs.append(newP)
                        p = newP

            pops.append(freqs)
            p = startingFreq

        for pop in pops:
            axs[ind].plot(pop)

        axs[ind].set_ylabel('Frequency of $p$')
        axs[ind].set_xlabel('Generations')
        axs[ind].set_ylim((0, 1))
        axs[ind].set_xlim((0, gens))
        axs[ind].set_title(f'2N = {twoN}')


    plt.show()
