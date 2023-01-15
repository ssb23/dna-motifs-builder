import numpy as np
from language import nucleotides
from language import converse
from basePenalties import BasePenalties


class BasePayloadPenalties(BasePenalties):
    def __init__(self, constraints, hyperparams):
        BasePenalties.__init__(self, constraints, hyperparams)
        self.payloads = set()
        self.motifSize = constraints.motifSize

        # Hairpin
        self.maxHairpin = constraints.maxHairpin

        # Homopolymer
        self.maxHom = constraints.maxHom
        self.maxStartHom = {}
        self.maxEndHom = {}

        # GC-content
        self.minGc = constraints.minGc
        self.maxGc = constraints.maxGc

        # Joint
        self.joints = set()
        self.startJoints = set()
        self.endJoints = set()
        self.keySize = constraints.keySize
        self.startJointHom = {'A':[], 'T':[], 'C':[], 'G':[]}
        self.endJointHom = {'A':[], 'T':[], 'C':[], 'G':[]}
        self.wholeJointHom = {'A':0, 'T':0, 'C':0, 'G':0}
        self.minJointsGcCount = -1
        self.maxJointsGcCount = -1

        # Payload
        self.payloadSize = constraints.payloadSize
        self.payloadNum = constraints.payloadNum
        self.startHomPayload = {'A':[], 'T':[], 'C':[], 'G':[]}
        self.endHomPayload = {'A':[], 'T':[], 'C':[], 'G':[]}

        # Hyperparameters
        self.homHyperparams = hyperparams.hom
        self.hairpinHyperparams = hyperparams.hairpin
        self.jointRepeatHyperparams = hyperparams.jointRepeat
        self.keyInPayloadHyperparams = hyperparams.keyInPayload
        self.keyGcContentHyperparams = hyperparams.keyGcContent
        self.motifGcContentHyperparams = hyperparams.motifGcContent

    ### Add Joints ###

    async def add_joints(self, joints):
        self.joints = joints
        await self.generate_joint_pre_stats()

    ### Pre-stats ###

    async def generate_joint_pre_stats(self):
        jointSize = int(self.keySize / 2)
        startMinGcCount = -1
        startMaxGcCount = -1
        endMinGcCount = -1
        endMaxGcCount = -1

        for joint in self.joints:
            cur = joint[0]
            curHom = 1
            curGcCount = 1 if cur in ['G', 'C'] else 0
            isEnd = True

            for i in range(1, self.keySize):
                nuc = joint[i]

                # Update end hom
                if cur != nuc:
                    if isEnd:
                        isEnd = False
                        self.endJointHom[cur].append(curHom)
                    cur = nuc
                    curHom = 0
                
                # Update end gc count
                if i == jointSize:
                    endMinGcCount = curGcCount if endMinGcCount == -1 else min(endMinGcCount, curGcCount)
                    endMaxGcCount =  max(endMaxGcCount, curGcCount)
                    curGcCount = 0
                
                curGcCount += 1 if nuc in ['G', 'C'] else 0

                curHom += 1

            # Update hom
            if isEnd:
                self.wholeJointHom[cur] += 1
            else:
                self.startJointHom[cur].append(curHom)

            # Update start gc count
            startMinGcCount = curGcCount if startMinGcCount == -1 else min(startMinGcCount, curGcCount)
            startMaxGcCount =  max(startMaxGcCount, curGcCount)
        # Update overall gc count
        self.minJointsGcCount = startMinGcCount + endMinGcCount
        self.maxJointsGcCount = startMaxGcCount + endMaxGcCount

    ### Get Penalties ###

    def get_all_penalties(self, motif, nucleotide):
        penalties = self.get_homopolymer_penalty(motif, nucleotide)
        penalties += self.get_hairpin_penalty(motif, nucleotide)
        penalties += self.get_gc_penalty(motif, nucleotide)
        return penalties

    def get_homopolymer_penalty(self, motif, nucleotide):
        curMotif = motif + nucleotide
        return self.homopolymer_penalty(curMotif)

    def get_hairpin_penalty(self, motif, nucleotide):
        curMotif = motif + nucleotide
        return self.hairpin_penalty(curMotif)

    def get_gc_penalty(self, motif, nucleotide):
        curMotif = motif + nucleotide
        return self.motif_GC_content_penalty(curMotif)

    ### Add Payloads ###
    
    async def add_payload(self, newPayload):
        self.payloads.add(newPayload)
        await self.add_hom_stats(newPayload)
    
    async def add_payloads(self, newPayloads):
        for newPayload in newPayloads:
            await self.add_payload(newPayload)

    ### Add pre-stats ###

    async def add_hom_stats(self, payload):
        isStart = True
        cur = payload[0]
        homCount = 1
        for i in range(1, self.payloadSize):
            nuc = payload[i]

            # Update homopolymer stats
            if cur != nuc:
                if isStart:
                    isStart = False
                    self.startHomPayload[cur].append(homCount)
                homCount = 0
                cur = nuc
            homCount += 1

        # Update homopolymer stats
        if isStart:
            self.startHomPayload[cur].append(homCount)
        self.endHomPayload[cur].append(homCount)
       
    ### Homopolymer penalty ###

    def homopolymer_penalty(self, curPayload):
        homStats = self.homopolymer_stats(curPayload)
        return homStats

    def homopolymer_stats(self, curPayload):
        curHoms = []
        curHom = 1
        cur = curPayload[-1]
        for i in range(len(curPayload) - 2, -1, -1):
            nuc = curPayload[i]
            if cur != nuc:
                break
            curHom += 1
        newNuc = curPayload[-1]

        if not self.startJointHom[newNuc] and not self.endJointHom[newNuc] and self.wholeJointHom[newNuc] == 0:
            return self.homHyperparams.hyperparameter**(curHom / self.maxHom)

        # start and end joints
        if curHom == len(curPayload) and len(curPayload) == self.payloadSize:
            if not self.startJointHom[newNuc] and self.wholeJointHom[newNuc] == 0:
                for endHom in self.endJointHom[newNuc]:
                    curHoms.append(self.homHyperparams.hyperparameter**((curHom + endHom)/self.maxHom))
            elif not self.endJointHom[newNuc] and self.wholeJointHom[newNuc] == 0:
                for startHom in self.startJointHom[newNuc]:
                    curHoms.append(self.homHyperparams.hyperparameter**((curHom + startHom)/self.maxHom))
            else:
                for startHom in self.startJointHom[newNuc]:
                    for endHom in self.endJointHom[newNuc]:
                        curHoms.append(self.homHyperparams.hyperparameter**((curHom + startHom + endHom)/self.maxHom))

            if self.wholeJointHom[newNuc] != 0:
                # Worst case: whole motif is homopolymer
                curHoms.append(self.homHyperparams.hyperparameter**(self.motifSize))

        # start joints
        elif curHom == len(curPayload):
            if not self.startJointHom[newNuc] and self.wholeJointHom[newNuc] == 0:
                curHoms.append(self.homHyperparams.hyperparameter**(curHom/self.maxHom))
            for startHom in self.startJointHom[newNuc]:
                curHoms.append(self.homHyperparams.hyperparameter**((curHom + startHom)/ self.maxHom))

            if self.wholeJointHom[newNuc] != 0:
                curHoms.append(self.homHyperparams.hyperparameter**((curHom + self.keySize)/ self.maxHom))
                
                for endHom in self.endHomPayload[newNuc]:
                    for i in range(self.wholeJointHom[newNuc]):
                        curHoms.append(self.homHyperparams.hyperparameter**((curHom + endHom + self.keySize)/ self.maxHom))
        
        # end joints
        elif len(curPayload) == self.payloadSize:
            if not self.endJointHom[newNuc] and self.wholeJointHom[newNuc] == 0:
                curHoms.append(self.homHyperparams.hyperparameter**(curHom/self.maxHom))
            for endHom in self.endJointHom[newNuc]:
                curHoms.append(self.homHyperparams.hyperparameter**((curHom + endHom)/ self.maxHom))

            if self.wholeJointHom[newNuc] != 0:
                curStart = newNuc
                curStartHom = 0
                for i in range(len(curPayload)):
                    nuc = curPayload[i]
                    if curStart != nuc:
                        break
                    curStartHom += 1
                    
                curHoms.append(self.homHyperparams.hyperparameter**((curHom + curStartHom + self.keySize)/ self.maxHom))
                for startHom in self.startHomPayload[newNuc]:
                    for i in range(self.wholeJointHom[newNuc]):
                        curHoms.append(self.homHyperparams.hyperparameter**((curHom + startHom + self.keySize)/ self.maxHom))
        
        else:
            curHoms.append(self.homHyperparams.hyperparameter**(curHom/ self.maxHom))
        return np.sum(np.array(curHoms))

    ### Motif GC Content Penalty ###

    def motif_GC_content_penalty(self, curPayload):
        motifGcContent = self.motif_GC_content_stats(curPayload)
        return motifGcContent

    def motif_GC_content_stats(self, curPayload):
        gcCount = 0
        for i in range(len(curPayload)):
            nuc = curPayload[i]
            gcCount += 1 if nuc in ['G', 'C'] else 0
        
        curMotifSize = len(curPayload) + self.keySize
        weight = curMotifSize / self.motifSize
        minGcContent = (100 * (gcCount + self.minJointsGcCount)) / curMotifSize
        maxGcContent = (100 * (gcCount + self.maxJointsGcCount)) / curMotifSize
        score = 0
        score = max(score, weight * (self.minGc - minGcContent))
        score = max(score, weight * (maxGcContent - self.maxGc))
        return score

    ### Hairpin Penalty ###

    def hairpin_penalty(self, curPayload):
        hairpinStats = self.forward_hairpin_counter(curPayload, self.payloads, self.joints) + self.backward_hairpin_counter(curPayload, self.payloads, self.joints)
        return hairpinStats


if __name__ == '__main__':
    basePenalties = BasePenalties()
