import numpy as np
from language import nucleotides
from language import converse
from basePenalties import BasePenalties


class BaseKeyPenalties(BasePenalties):
    def __init__(self, constraints, hyperparams):
        BasePenalties.__init__(self, constraints, hyperparams)
        self.keys = set()
        self.keySize = constraints.keySize

        # Hairpin
        self.maxHairpin = constraints.maxHairpin

        # Homopolymer
        self.maxHom = constraints.maxHom
        self.curHom = 0

        # GC-content
        self.minGc = constraints.minGc
        self.maxGc = constraints.maxGc
        self.curGcCount = 0

        # Joints
        self.startJoints = set()
        self.endJoints = set()

        # Key
        self.keySize = constraints.keySize
        self.startMinGcCountMotif = -1
        self.endMinGcCountMotif = -1
        self.startMaxGcCountMotif = -1
        self.endMaxGcCountMotif = -1
        self.curKey = ''

        # Hyperparameters
        self.homHyperparams = hyperparams.hom
        self.hairpinHyperparams = hyperparams.hairpin
        self.keyGcContentHyperparams = hyperparams.keyGcContent
        self.similarityHyperparams = hyperparams.similarity
        self.uniqueJointsHyperparams = hyperparams.uniqueJoints

    ### Get Penalties ###

    def get_all_penalties(self, key, nucleotide):
        penlties = self.get_homopolymer_penalty(key, nucleotide)
        penlties += self.get_hairpin_penalty(key, nucleotide)
        penlties += self.get_gc_penalty(key, nucleotide)
        penlties += self.get_similarity_penalty(key, nucleotide)
        penlties += self.get_unique_joints_penalty(key, nucleotide)
        return penlties

    def get_homopolymer_penalty(self, key, nucleotide):
        curKey = key + nucleotide
        return self.homopolymer_penalty(curKey)

    def get_hairpin_penalty(self, key, nucleotide):
        curKey = key + nucleotide
        return self.hairpin_penalty(curKey)

    def get_gc_penalty(self, key, nucleotide):
        curKey = key + nucleotide
        return self.key_GC_content_penalty(curKey)
    
    def get_similarity_penalty(self, key, nucleotide):
        curKey = key + nucleotide
        return self.similarity_penalty(curKey)

    def get_unique_joints_penalty(self, key, nucleotide):
        curKey = key + nucleotide
        return self.unique_joints_penalty(curKey)
    
    ### Add Base to current Key ###
    
    async def add_base(self, newBase):
        self.curKey += newBase
        await self.update_homopolymer_stats(newBase)
        await self.update_gc_count_stats(newBase)

    async def update_homopolymer_stats(self, newBase):
        if len(self.curKey) == 1:
            self.curHom = 1
            return
        if self.curKey[len(self.curKey) - 2] != newBase:
            self.curHom = 0
        self.curHom += 1 

    async def update_gc_count_stats(self, newBase):
        self.curGcCount += 1 if newBase in ['G', 'C'] else 0

    ### Start New Key ###

    def start_new_key(self):
        self.curKey = ''
        self.curGcCount = 0
        self.curHom = 0

    ### Add Keys ###

    async def add_key(self, newKey):
        if len(newKey) != self.keySize:
            return False
        self.keys.add(newKey)
        await self.add_joints(newKey)
        await self.add_motif_gc_info(newKey)
    
    async def add_keys(self, newKeys):
        for newKey in newKeys:
            await self.add_key(newKey)
        import time

    async def add_joints(self, newKey):
        jointSize = int(self.keySize / 2)
        self.startJoints.add(newKey[:jointSize])
        self.endJoints.add(newKey[jointSize:])

    async def add_motif_gc_info(self, k):
        jointSize = int(self.keySize / 2)
        startMinGcCount = -1
        startMaxGcCount = -1
        endMinGcCount = -1
        endMaxGcCount = -1

        cur = k[0]
        curGcCount = 1 if cur in ['G', 'C'] else 0

        for i in range(1, len(k)):
            nuc = k[i]

            # Update start gc count
            if i == jointSize:
                startMinGcCount = curGcCount if startMinGcCount == -1 else min(startMinGcCount, curGcCount)
                startMaxGcCount =  max(startMaxGcCount, curGcCount)
                curGcCount = 0
            
            # motif gc is complementary
            curGcCount += 1 if nuc in ['G', 'C'] else 0

        # Update end gc count
        endMinGcCount = curGcCount if endMinGcCount == -1 else min(endMinGcCount, curGcCount)
        endMaxGcCount =  max(endMaxGcCount, curGcCount)

        # Update overall gc count
        self.startMinGcCountMotif = startMinGcCount if self.startMinGcCountMotif == -1 else min(self.startMinGcCountMotif, startMinGcCount)
        self.endMinGcCountMotif = endMinGcCount if self.endMinGcCountMotif == -1 else min(self.endMinGcCountMotif, endMinGcCount)
        self.startMaxGcCountMotif = max(self.startMaxGcCountMotif, startMaxGcCount)
        self.endMaxGcCountMotif = max(self.endMaxGcCountMotif, endMaxGcCount)

    ### Homopolymer penalty ###

    def homopolymer_penalty(self, curKey):
        homSize = self.homopolymer_stats(curKey)
        if homSize == 0:
            return 0
        return self.homHyperparams.hyperparameter**(homSize / self.maxHom)

    def homopolymer_stats(self, curKey):
        if len(curKey) <= 1:
            return len(curKey)
        return self.curHom + 1 if curKey[len(curKey) - 1] == curKey[len(curKey) - 2] else 1

    ### Key GC Content Penalty ###

    def key_GC_content_penalty(self, curKey):
        keyGcContent = self.key_GC_content_stats(curKey)
        return keyGcContent

    def key_GC_content_stats(self, curKey):
        score = self.key_GC_content_within_key(curKey) + self.key_GC_content_within_motif(curKey)
        return score

    def key_GC_content_within_key(self, curKey):
        gcCount = self.curGcCount + 1 if curKey[len(curKey) - 1] in ['G', 'C'] else self.curGcCount
        gcContent = (100 * gcCount) / len(curKey)
        penalty = 0
        penalty = max(penalty, gcContent - self.maxGc)
        penalty = max(penalty, self.minGc - gcContent)
        return penalty * (len(curKey) / self.keySize)

    def key_GC_content_within_motif(self, curKey):
        keyGcCount = self.curGcCount + 1 if curKey[len(curKey) - 1] in ['G', 'C'] else self.curGcCount
        jointSize = int(self.keySize / 2)

        # motif gc is complementary
        if len(curKey) <= jointSize:
            startGcCountMotif = keyGcCount
            curMotifSize = len(curKey) + self.payloadSize + jointSize
            weight = curMotifSize / self.motifSize
            minGcContent = (100 * (startGcCountMotif + max(self.endMinGcCountMotif, 0))) / curMotifSize
            maxGcContent = (100 * (startGcCountMotif + max(self.endMaxGcCountMotif, 0))) / curMotifSize
        else:
            endGcCountMotif = sum([1 if curKey[x] in ['G', 'C'] else 0 for x in range(jointSize, len(curKey))])
            startGcCountMotif = keyGcCount - endGcCountMotif
            startMinGcCountMotif = startGcCountMotif if self.startMinGcCountMotif == -1 else min(self.startMinGcCountMotif, startGcCountMotif)
            startMaxGcCountMotif = max(self.startMaxGcCountMotif, startGcCountMotif)
            curMotifSize = len(curKey) + self.payloadSize
            weight = curMotifSize / self.motifSize
            minGcContent = (100 * (endGcCountMotif + startMinGcCountMotif)) / curMotifSize
            maxGcContent = (100 * (endGcCountMotif + startMaxGcCountMotif)) / curMotifSize
        score = 0
        score = max(score, weight * (self.minGc - minGcContent))
        score = max(score, weight * (maxGcContent - self.maxGc))
        return score

    ### Hairpin Penalty ###

    def hairpin_penalty(self, curKey):
        # since (complementary hairpin) == hairpin
        hairpinStats = self.forward_hairpin_counter(curKey, self.keys, set(), isKey=True) + self.backward_hairpin_counter(curKey, self.keys, set(), isKey=True)
        return hairpinStats

    ### Unique Joints ###
    
    def unique_joints_penalty(self, curKey):
        uniqueJointsStats = self.unique_joints_stats(curKey)
        return uniqueJointsStats

    def unique_joints_stats(self, curKey):
        jointSize = int(self.keySize / 2)
        numPossibleUniqueJoints = 4**jointSize
        if len(self.startJoints) == numPossibleUniqueJoints and len(self.endJoints) == numPossibleUniqueJoints:
            return 0
        maxSimilarity = 0
        weight = 0

        if len(curKey) <= jointSize:
            for startJoint in self.startJoints:
                if curKey == startJoint[:len(curKey)]:
                    maxSimilarity = max(maxSimilarity, len(curKey))
            weight = 100 *(numPossibleUniqueJoints - len(self.startJoints)) / numPossibleUniqueJoints
        else:
            for endJoint in self.endJoints:
                if curKey[jointSize:] == endJoint[:len(curKey) - jointSize]:
                    maxSimilarity = max(maxSimilarity, len(curKey) - jointSize)
            weight = 100 * (numPossibleUniqueJoints - len(self.endJoints)) / numPossibleUniqueJoints
        return weight * self.uniqueJointsHyperparams.hyperparameter**(maxSimilarity / jointSize)

    ### Similarity ###

    def similarity_penalty(self, curKey):
        similarityStats = self.similarity_stats(curKey)
        return similarityStats

    def similarity_stats(self, curKey):
        maxSimilarity = 0
        
        for k in self.keys:
            if curKey == k[:len(curKey)]:
                maxSimilarity = max(maxSimilarity, len(curKey))
        return self.similarityHyperparams.hyperparameter**(maxSimilarity / self.keySize)
    

if __name__ == '__main__':
    basePenalties = BaseKeyPenalties()
