from basePayloadPenalties import BasePayloadPenalties
from constraints import Constraints
from validation import Validate
from hyperparameters import Hyperparameters
from language import nucleotides

import numpy as np


class MotifBuilder:
    def __init__(self, constraints, hyperparameters, weights):
        self.keySize = constraints.keySize
        self.constraints = constraints

        self.hyperparameters = hyperparameters
        self.weights = weights

        self.payloadNum = constraints.payloadNum
        self.payloadSize = constraints.payloadSize
       # self.motifSize = constraints.motifSize

        self.penalties = BasePayloadPenalties(constraints, hyperparameters)
        
    async def add_joints(self, joints):
        await self.penalties.add_joints(joints)

    ### Build Payload ###

    async def buildPayload(self, withConstraints):
        payload = ''

        for _ in range(self.payloadSize):
            
            constraintPenalties = {}
            # 'hom', 'hairpin', 'motifGcContent'
            for constraint in withConstraints:
                constraintPenalties[constraint] = []

            for n in nucleotides:
                for constraint in constraintPenalties:
                    if constraint == 'hom':
                        homPenalty = self.penalties.get_homopolymer_penalty(payload, n)
                        constraintPenalties[constraint].append(homPenalty)
                    elif constraint == 'hairpin':
                        hairpinPenalty = self.penalties.get_hairpin_penalty(payload, n)
                        constraintPenalties[constraint].append(hairpinPenalty)
                    elif constraint == 'motifGcContent':
                        motifGcPenalty = self.penalties.get_gc_penalty(payload, n)
                        constraintPenalties[constraint].append(motifGcPenalty)
            
            for constraint in constraintPenalties:
                constraintPenalties[constraint] = np.array(constraintPenalties[constraint])
            
            if not withConstraints:
                p = [0.25, 0.25, 0.25, 0.25]
            else:
                lp = 0
                for constraint in constraintPenalties:
                    weight = self.weights[constraint]
                    lp -= 1.0 * weight * constraintPenalties[constraint]
                p = np.exp(lp)
                p /= p.sum()

            new_nucleotide = np.random.choice(nucleotides, p=p)

            payload += new_nucleotide

        await self.penalties.add_payload(payload)

        return payload
    
    async def buildAllPayloads(self, withConstraints):
        allPayloads = set()
        for _ in range(self.payloadNum):
            payload = await self.buildPayload(withConstraints)
            allPayloads.add(payload)
        return allPayloads


async def main():
    payloadSize = 8
    payloadNum = 5
    maxHom = 1
    maxHairpin = 2
    loopSize = 1
    minGc = 25
    maxGc = 60
    keySize = 2
    keyNum = 4
    
    constraints = Constraints(payloadSize, payloadNum, maxHom, maxHairpin, loopSize, minGc, maxGc, keySize, keyNum)

    import matplotlib.pyplot as plt

    y = []
    allKeysScore = 0
    numRounds = 1
    numBadSequences = 0
    withConstraints = {'hom', 'motifGcContent', 'hairpin'} #'similarity', 'uniqueJoints'}
    for i in range(numRounds):
        hyperparams = {'hom': 5, 'motifGcContent': 5, 'hairpin': 10}
        weights = {'hom': 1, 'motifGcContent': 1, 'hairpin': 1}
        internalHyperparameters = Hyperparameters(hyperparams)

        keys =  {'GT', 'GA', 'CT', 'TG'} # {'AT', 'GC'}
        joints = {'AC', 'TC', 'AG', 'CA'}

        motifBuilder = MotifBuilder(constraints, internalHyperparameters, weights)
        await motifBuilder.add_joints(joints)
        payloads = await motifBuilder.buildAllPayloads(withConstraints)
        
        print(payloads)

        payloadsValidation = Validate(constraints)
        payloadsValidation.add_keys(keys)
        await payloadsValidation.add_payloads(payloads)

        homScore = payloadsValidation.get_motifs_and_keys_homopolymer_score()
        gcScore = payloadsValidation.get_motifs_gc_score()
        hairpinScore = payloadsValidation.get_motifs_and_keys_hairpin_score()
        totalScore = homScore + gcScore + hairpinScore
        print('homScore: ', homScore)
        print('gcScore: ', gcScore)
        print('hairpinScore: ', hairpinScore)

        allKeysScore += totalScore
        numBadSequences += 0 if totalScore == 0 else 1
        y.append(totalScore)
    print('Total score of all keys sets is ', allKeysScore)
    print('Number of bad keys sets is ', numBadSequences)
    x = range(len(y))
    plt.plot(x,y)
   # plt.show()


if __name__ == '__main__':
    import asyncio
    asyncio.run(main())
