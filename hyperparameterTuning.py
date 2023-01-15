from constraints import Constraints
from validation import Validate
from hyperparameters import Hyperparameters
from motifBuilder import MotifBuilder

import numpy as np
import scipy.optimize
from time import time


class HyperparameterTuning:
    def __init__(self, constraints, constraints_to_tune, number_of_tries=1000):
        self.constraints = constraints
        self.initial_guess = []
        self.final_guess = []
        self.number_of_tries = number_of_tries
        self.withConstraints = {}

        self.y = []

        self.timeTaken = 0  

    def set_number_of_tries(self, num_tries):
        self.number_of_tries = num_tries

    def set_initial_hyperparameters(self, constraints_to_tune):
        self.initial_guess = []
        for constraint in constraints_to_tune:
            # weight
            self.initial_guess.append(1)
            # hyperparameter
            self.initial_guess.append(5)
        self.initial_guess = np.array(self.initial_guess)

    def get_time_taken(self):
        return self.timeTaken

    def get_initial_hyperparameters(self):
        return self.initial_guess
    
    def get_hyperparameters(self):
        return self.final_guess

    def run_hyperparameter_tuning(self, constraints_to_tune):
        self.set_initial_hyperparameters(constraints_to_tune)
        self.constraints_to_tune = constraints_to_tune
        t0 = time()
        # Bounds: want all hyperparameters to be > 1 
        maximum_a_posteriori = scipy.optimize.minimize(
            lambda x: -self.log_posterior(x), self.initial_guess, \
             bounds=[(1,None)]*len(self.initial_guess))
        self.timeTaken = time() - t0

        self.final_guess = maximum_a_posteriori.x
        return maximum_a_posteriori.x

    def log_prior(self, hyperparameter_array):
        # Add any constraints on hyperparameters here
        return 0

    def quality_check(self, set_of_motifs):
        # 0 is good, score < 0 if bad
        quality = Validate(set_of_motifs, self.constraints)
        quality.calculate_scores_of_constraints(constraints_to_tune)
    
        return quality.get_total_score()

    def f(self, hyperparameter_array):
        hyperparams = {}
        weights = {}
        i = 0
        for constraint in self.constraints_to_tune:
            weights[constraint] = hyperparameter_array[i]
            hyperparams[constraint] = hyperparameter_array[i + 1]
            i += 1

        hyperparameters = Hyperparameters(hyperparams)

        motifBuilder = MotifBuilder(self.constraints, hyperparameters, weights)
        motifBuilder.buildAllMotifs(self.constraints_to_tune)
        return motifBuilder.getAllMotifs()

    def log_likelihood(self, hyperparameter_array):
        quality = 0
       
        for k in range(self.number_of_tries):
            set_of_motifs = self.f(hyperparameter_array) # always different results
            q = self.quality_check(set_of_motifs)
            quality += q
            self.y.append(q)
        quality /= self.number_of_tries

        return quality
    
    def log_posterior(self, hyperparameter_array):
        return self.log_prior(hyperparameter_array) + self.log_likelihood(hyperparameter_array)


if __name__ == '__main__':
    motifSize = 20
    motifNum = 10
    maxHom = 2
    maxHairpin = 2
    loopSize = 1
    minGc = 25
    maxGc = 60
    keySize = 2
    
    constraints = Constraints(motifSize, motifNum, maxHom, maxHairpin, loopSize, minGc, maxGc, keySize)
    
    import matplotlib.pyplot as plt
    
    constraints_to_tune = {'hom'}
    number_of_tries = 1

    h = HyperparameterTuning(constraints, number_of_tries)
    h.run_hyperparameter_tuning(constraints_to_tune)

    print(h.final_guess)
    print(h.get_time_taken())
    x = range(len(h.y))
    plt.plot(x,h.y)
    plt.show()
