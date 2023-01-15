from baseKeyPenalties import BaseKeyPenalties
from constraints import Constraints
from validation import Validate
from hyperparameters import Hyperparameters
from language import nucleotides

import numpy as np


class KeyBuilder:
    def __init__(self, constraints, hyperparameters, weights):
        self.constraints = constraints
        
        self.hyperparameters = hyperparameters
        self.weights = weights

        self.keyNum = constraints.keyNum
        self.keySize = constraints.keySize

        self.penalties = BaseKeyPenalties(constraints, hyperparameters)

    ### Build Keys ###

    async def buildKey(self, withConstraints):
        key = ''

        self.penalties.start_new_key()
        for _ in range(self.keySize):
            
            constraintPenalties = {}
            # 'hom', 'hairpin', 'keyGcContent'
            for constraint in withConstraints:
                constraintPenalties[constraint] = []

            for n in nucleotides:

                for constraint in constraintPenalties:
                    if constraint == 'hom':
                            homPenalty = self.penalties.get_homopolymer_penalty(key, n)
                            constraintPenalties[constraint].append(homPenalty)
                    if constraint == 'hairpin':
                            hairpinPenalty = self.penalties.get_hairpin_penalty(key, n)
                            constraintPenalties[constraint].append(hairpinPenalty)
                    if constraint == 'keyGcContent':
                            keyGcPenalty = self.penalties.get_gc_penalty(key, n)
                            constraintPenalties[constraint].append(keyGcPenalty)
                        # case 'similarity':
                        #     similarity = self.get_similarity_penalty(key, n)
                        #     constraintPenalties[constraint].append(similarity)
                        # case 'uniqueJoints':
                        #     uniqueJoints = self.get_unique_joints_penalty(key, n)
                        #     constraintPenalties[constraint].append(uniqueJoints)
            
                

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
        
            key += new_nucleotide
            await self.penalties.add_base(new_nucleotide)

        await self.penalties.add_key(key)
        return key
    
    async def buildAllKeys(self, withConstraints):
        allKeys = set()
        for _ in range(self.keyNum):
            key = await self.buildKey(withConstraints)
            allKeys.add(key)
        return allKeys


async def main():
    payloadSize = 8
    payloadNum = 5
    maxHom = 1
    maxHairpin = 2
    loopSize = 1
    minGc = 25
    maxGc = 60
    keySize = 2
    keyNum = 5
    
    constraints = Constraints(payloadSize, payloadNum, maxHom, maxHairpin, loopSize, minGc, maxGc, keySize, keyNum)

    import matplotlib.pyplot as plt

    y = []
    allKeysScore = 0
    numRounds = 1
    numBadSequences = 0
    withConstraints = {'hom', 'keyGcContent', 'hairpin'} #'similarity', 'uniqueJoints'}
    for i in range(numRounds):
        hyperparams = {'hom': 5, 'keyGcContent': 5, 'hairpin': 5, 'similarity': 5, 'uniqueJoints': 5}
        weights = {'hom': 1, 'keyGcContent': 1, 'hairpin': 1, 'similarity': 1, 'uniqueJoints': 1}
        internalHyperparameters = Hyperparameters(hyperparams)

        keyBuilder = KeyBuilder(constraints, internalHyperparameters, weights)
        keys = await keyBuilder.buildAllKeys(withConstraints)
        
        print(keys)

        keysValidation = Validate(constraints)
        keysValidation.add_keys(keys)
        homScore = keysValidation.get_keys_homopolymer_score()
        gcScore = keysValidation.get_keys_gc_score()
        hairpinScore = keysValidation.get_keys_hairpin_score()
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
    
#     payloadSize = 8
#     payloadNum = 5
#     maxHom = 1
#     maxHairpin = 2
#     loopSize = 1
#     minGc = 25
#     maxGc = 60
#     keySize = 2
#     keyNum = 5
    
#     constraints = Constraints(payloadSize, payloadNum, maxHom, maxHairpin, loopSize, minGc, maxGc, keySize, keyNum)

#     import matplotlib.pyplot as plt

#     y = []
#     allKeysScore = 0
#     numRounds = 1
#     numBadSequences = 0
#     withConstraints = {'keyHom', 'keyGcContent', 'keyHairpin'} #'similarity', 'uniqueJoints'}
#     for i in range(numRounds):
#         hyperparams = {'hom': 5, 'keyGcContent': 5, 'hairpin': 5, 'similarity': 5, 'uniqueJoints': 5}
#         weights = {'hom': 1, 'keyGcContent': 1, 'hairpin': 1, 'similarity': 1, 'uniqueJoints': 1}
#         internalHyperparameters = Hyperparameters(hyperparams)

#         keyBuilder = KeyBuilder(constraints, internalHyperparameters, weights)
#         keyBuilder.buildAllKeys(withConstraints)
#         keys = 

#         print(keys)
#         print('MOSMDOSAMDOSMDOSAMDOMSDOMSAODMOS', keys)


#         keysValidation = Validate(constraints)
#         keysValidation.add_keys(keys)
#         homScore = keysValidation.get_keys_homopolymer_score()
#         gcScore = keysValidation.get_keys_gc_score()
#         hairpinScore = keysValidation.get_keys_hairpin_score()
#         totalScore = homScore + gcScore + hairpinScore
#         print('homScore: ', homScore)
#         print('gcScore: ', gcScore)
#         print('hairpinScore: ', hairpinScore)

#         allKeysScore += totalScore
#         numBadSequences += 0 if totalScore == 0 else 1
#         y.append(totalScore)
#     print('Total score of all keys sets is ', allKeysScore)
#     print('Number of bad keys sets is ', numBadSequences)
#     x = range(len(y))
#     plt.plot(x,y)
#    # plt.show()

