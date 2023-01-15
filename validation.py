from language import nucleotides
from language import converse
import numpy as np
from basePenalties import BasePenalties

class Validate (BasePenalties):
    # 0 is good, score < 0 if bad
    def __init__(self, constraints):
        BasePenalties.__init__(self, constraints)

        # Keys
        self.keys = set()

        # Payloads
        self.payloads = set()

        # Joints
        self.joints = set()
        self.startJointHom = {'A':[], 'T':[], 'C':[], 'G':[]}
        self.endJointHom = {'A':[], 'T':[], 'C':[], 'G':[]}
        self.wholeJointHom = {'A':0, 'T':0, 'C':0, 'G':0}

        self.wholeHomPayload = {'A':0, 'T':0, 'C':0, 'G':0}
        self.startHomPayload = {'A':[], 'T':[], 'C':[], 'G':[]}
        self.endHomPayload = {'A':[], 'T':[], 'C':[], 'G':[]}

        self.startJointsGcCount = []
        self.endJointsGcCount = []
        self.payloadsGcCount = []

        # Motifs
        self.motifSize = constraints.motifSize
        self.maxPossibleMotifNum = constraints.maxPossibleMotifNum

        # Hairpin
        self.maxHairpin = constraints.maxHairpin
        self.loopSizeMin = constraints.loopSizeMin
        self.loopSizeMax = constraints.loopSizeMax

        # Homopolymer
        self.maxHom = constraints.maxHom
        self.homOutBoundary = []

        # GC-content
        self.minGc = constraints.minGc
        self.maxGc = constraints.maxGc
        self.gcContentsOutBoundary = []

        self.payloadsGcCount = []

        # Key
        self.startKeys = set()
        self.endKeys = set()
        self.keySize = constraints.keySize
        # GC-content key
        self.startJointsGcCount = []
        self.endJointsGcCount = []
        self.keyGcContentsOutBoundary = []

        # Penalties
        self.motifHomPenalty = 0
        self.motifGcPenalty = 0
        self.motifHairpinPenalty = 0
        self.keyHomPenalty = 0
        self.keyGcPenalty = 0
        self.keyHairpinPenalty = 0

        self.jointRepeatPenalty = 0
        self.keyInPayloadPenalty = 0
        
        # Total score 0 is no penalty and negatif is penalty
        self.score = 0

    ### Add Keys and Payloads ###

    async def add_keys_and_payloads(self, keys, payloads):
        self.keys = keys
        self.payloads = payloads
        # Generate pre-stats
        await self.generate_motifs_stats()

    async def add_payloads(self, payloads):
        self.payloads = payloads
        # Generate pre-stats
        await self.generate_motifs_stats()

    def add_keys(self, keys):
        self.keys = keys

    ### Calculate Scores ###

    def get_all_motifs_and_keys_scores_total(self):
        score = self.get_motifs_and_keys_homopolymer_score()
        score += self.get_motifs_and_keys_hairpin_score()
        score += self.get_motifs_gc_score()
        score += self.get_keys_gc_score()
        return score
    
    def get_all_keys_scores_total(self):
        score = self.get_keys_homopolymer_score()
        score += self.get_keys_hairpin_score()
        score += self.get_keys_gc_score()
        return score
    
    def get_total_scores_of_constraints(self, withConstraints):
        score = 0
        for constraint in withConstraints:
            if constraint == 'hom':
                    score += self.get_motifs_and_keys_homopolymer_score()
            if constraint == 'keyHom':
                    score += self.get_keys_homopolymer_score()
            if constraint == 'hairpin':
                    score += self.get_motifs_and_keys_hairpin_score()
            if constraint == 'keyHairpin':
                    score += self.get_keys_hairpin_score()
            if constraint == 'motifGcContent':
                    score += self.get_motifs_gc_score()
            if constraint == 'keyGcContent':
                    score += self.get_keys_gc_score()
        return score

    def get_motifs_and_keys_homopolymer_score(self):
        return -self.get_motif_and_key_homopolymer_penalty()

    def get_keys_homopolymer_score(self):
        return -self.get_key_homopolymer_penalty()
    
    def get_motifs_and_keys_hairpin_score(self):
        return -self.get_motif_and_key_hairpin_penalty()

    def get_keys_hairpin_score(self):
        return -self.get_key_hairpin_penalty()

    def get_motifs_gc_score(self):
        return -self.get_motif_gc_penalty()

    def get_keys_gc_score(self):
        return -self.get_key_gc_penalty()

    ##### Key Penalties Only #####

    ### Generate Key Homopoylmer Score (no need to use if check motif homopolymers) ###

    def get_key_homopolymer_penalty(self):
        homs = []
        for k in self.keys:
            curBase = k[0]
            hom = 1
            for i in range(1, self.keySize):
                nuc = k[i]
                if nuc != curBase:
                    if hom - self.maxHom > 0:
                        homs.append(hom - self.maxHom)
                    curBase = nuc
                    hom = 0
                hom += 1
            if hom - self.maxHom > 0:
                homs.append(hom - self.maxHom)
        return np.sum(np.array(homs))

    ### Generate Key GC Score ###

    def get_key_gc_penalty(self):
        gcScores = []
        for k in self.keys:
            gcCount = 0
            for i in range(self.keySize):
                nuc = k[i]
                gcCount += 1 if nuc in ['G', 'C'] else 0

            gcContent = (100 * gcCount) / self.keySize
            score = 0
            score = max(score, self.minGc - gcContent)
            score = max(score, gcContent - self.maxGc)
            if score > 0:
                gcScores.append(score)

        return np.sum(np.array(gcScores))

    ### Generate Key Hairpin Score (no need to use if check motif hairpins) ###

    def get_key_hairpin_penalty(self):
        hairpinCount = self.key_hairpin_count()
        return hairpinCount

    def key_hairpin_count(self):
        hairpinCount = 0
        for loopSize in range(self.loopSizeMin, self.loopSizeMax + 1):
            if self.maxHairpin * 2 + loopSize >= self.keySize:
                continue
            for key in self.keys:
                for stem1Start in range(self.keySize - self.maxHairpin * 2 - loopSize + 1):
                    hairpinCount += self.forward_hairpin_counter_at_startPos(key, self.keys, set(), stem1Start, [], isKey=True, loopSizeMin=loopSize, loopSizeMax=loopSize)
        return hairpinCount
    
    ##### Motif Penalties #####
    
    ### Generate Motif Homopolymer Score ###

    def get_motif_and_key_homopolymer_penalty(self):
        homOutBoundary = self.add_homopolymer_over_motifs_stats()
        if len(homOutBoundary) == 0:
            homStats = 0
        else:
            homStats = np.sum(np.array(homOutBoundary))
        return homStats

    def add_homopolymer_over_motifs_stats(self):
        jointSize = int(self.keySize / 2)
        
        for nuc in nucleotides:
            # worst case: whole motif is homopolymer
            if self.wholeHomPayload[nuc] > 0 and self.wholeJointHom[nuc] > 0:
                return [self.motifSize * self.maxPossibleMotifNum]
            # start and end joints
            if self.wholeHomPayload[nuc] > 0:
                if len(self.startJointHom[nuc]) == 0 and len(self.endJointHom[nuc]) == 0:
                    if self.payloadSize - self.maxHom > 0:
                        self.homOutBoundary.append(self.payloadSize - self.maxHom)
                elif len(self.startJointHom[nuc]) == 0:
                    for hom in self.endJointHom[nuc]:
                        if self.payloadSize + hom - self.maxHom > 0:
                            self.homOutBoundary.append(self.payloadSize + hom - self.maxHom)
                elif len(self.endJointHom[nuc]) == 0:
                    for hom in self.startJointHom[nuc]:
                        if self.payloadSize + hom - self.maxHom > 0:
                            self.homOutBoundary.append(self.payloadSize + hom - self.maxHom)
                else:
                    for hom1 in self.startHomPayload[nuc]:
                        for hom2 in self.endHomPayload[nuc]:
                            if self.payloadSize + hom1 + hom2 - self.maxHom > 0:
                                self.homOutBoundary.append(self.payloadSize + hom1 + hom2 - self.maxHom)
            # start and end payloads
            if self.wholeJointHom[nuc] > 0:
                if len(self.startHomPayload[nuc]) == 0 and len(self.endHomPayload[nuc]) == 0:
                    if self.keySize - self.maxHom > 0:
                        self.homOutBoundary.append(self.keySize - self.maxHom)
                elif len(self.startHomPayload[nuc]) == 0:
                    for hom in self.endHomPayload[nuc]:
                        if self.keySize + hom - self.maxHom > 0:
                            self.homOutBoundary.append(self.keySize + hom - self.maxHom)
                elif len(self.endHomPayload[nuc]) == 0:
                    for hom in self.startHomPayload[nuc]:
                        if self.keySize + hom - self.maxHom > 0:
                            self.homOutBoundary.append(self.keySize + hom - self.maxHom)
                else:
                    for hom1 in self.startHomPayload[nuc]:
                        for hom2 in self.endHomPayload[nuc]:
                            if self.keySize + hom1 + hom2 - self.maxHom > 0:
                                self.homOutBoundary.append(self.keySize + hom1 + hom2 - self.maxHom)
            # end payload with end joint
            for homPayload in self.endHomPayload[nuc]:
                for homJoint in self.endJointHom[nuc]:
                    if homPayload + homJoint - self.maxHom > 0:
                        self.homOutBoundary.append(homPayload + homJoint - self.maxHom)
            # start payload with start joint
            for homPayload in self.startHomPayload[nuc]:
                for homJoint in self.startJointHom[nuc]:
                    if homPayload + homJoint - self.maxHom > 0:
                        self.homOutBoundary.append(homPayload + homJoint - self.maxHom)
        homOutBoundary = self.homOutBoundary
        return homOutBoundary

    ### Generate Motif GC Score ###

    def get_motif_gc_penalty(self):
        gcScores = []
        for payloadGcCount in self.payloadsGcCount:
            for startJointGcCount in self.startJointsGcCount:
                for endJointGcCount in self.endJointsGcCount:
                    gcCount = payloadGcCount + startJointGcCount + endJointGcCount
                    gcContent = (100 * gcCount) / self.motifSize
                    score = 0
                    score = max(score, self.minGc - gcContent)
                    score = max(score, gcContent - self.maxGc)
                    if score > 0:
                        gcScores.append(score)
        return np.sum(np.array(gcScores))

    ### Generate Motif Hairpin Score ###

    def get_motif_and_key_hairpin_penalty(self):
        hairpinStats = self.hairpin_stats()
        return hairpinStats

    def hairpin_stats(self):  
        hairpinsCount = 0

        # forward hairpin
        
        for payload in self.payloads:
            for stem1Start in range(self.payloadSize):
                hairpinsCount += self.forward_hairpin_counter_at_startPos(payload, self.payloads, self.joints, stem1Start, [])

        for joint in self.joints:
            for stem1Start in range(self.keySize):
                hairpinsCount += self.forward_hairpin_counter_at_startPos(joint, self.joints, self.payloads, stem1Start, [], isKey=True)

        return hairpinsCount
    
    ##### Generate pre-stats #####
    
    async def generate_motifs_stats(self):
        jointSize = int(self.keySize / 2)

        # Joint stats
        # Joints = reverse complementary of keys
        for k in self.keys:
            cur = k[0]
            curHom = 1
            curGcCount = 1 if cur in ['G', 'C'] else 0
            isEnd = True
            joint = converse[k[0]]

            for i in range(1, self.keySize):
                nuc = k[i]
                joint = converse[nuc] + joint

                # Update start hom
                if cur != nuc:
                    if isEnd:
                        isEnd = False
                        self.startJointHom[converse[cur]].append(curHom)
                    elif curHom - self.maxHom > 0:
                        self.homOutBoundary.append(curHom - self.maxHom)
                    cur = nuc
                    curHom = 0
                
                # Update end gc count
                if i == jointSize:
                    self.endJointsGcCount.append(curGcCount)
                    curGcCount = 0
                
                curGcCount += 1 if nuc in ['G', 'C'] else 0

                curHom += 1
            
            self.joints.add(joint)

            # Update end hom
            if isEnd:
                self.wholeJointHom[converse[cur]] += 1
            else:
                self.endJointHom[converse[cur]].append(curHom)

            # Update start gc count
            self.startJointsGcCount.append(curGcCount)

        # Payload stats
        for payload in self.payloads:
            isStart = True
            cur = payload[0]
            homCount = 1
            gcCount = 1 if cur in ['G', 'C'] else 0
            for i in range(1, self.payloadSize):
                nuc = payload[i]

                # Update homopolymer stats
                if cur != nuc:
                    if isStart:
                        isStart = False
                        self.startHomPayload[cur].append(homCount)
                    elif homCount - self.maxHom > 0:
                        self.homOutBoundary.append(homCount - self.maxHom)
                    
                    homCount = 0
                    cur = nuc
                homCount += 1

                # Update GC Count
                if nuc in ['G', 'C']:
                    gcCount += 1

            # Update homopolymer stats
            if isStart:
                self.wholeHomPayload[cur] += homCount
            else:
                self.endHomPayload[cur].append(homCount)

            # Update GC Count
            self.payloadsGcCount.append(gcCount)


if __name__ == '__main__':
    motifs = {'ATCGGCGCGC', 'CAGTGATACGATCG'}
    validate = Validate(constraints, keys, payloads)