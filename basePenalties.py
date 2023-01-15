import numpy as np
from language import nucleotides
from language import converse


class BasePenalties:
    def __init__(self, constraints, hyperparams=None, joints=set(), payloads=set()):
        self.motifSize = constraints.motifSize

        # Hairpin
        self.maxHairpin = constraints.maxHairpin
        self.loopSizeMin = constraints.loopSizeMin
        self.loopSizeMax = constraints.loopSizeMax

        # Joints = reverse complementary of keys
        self.joints = joints

        # Key
        self.keySize = constraints.keySize

        # Payload
        self.payloads = payloads
        self.payloadSize = constraints.payloadSize

        # Penalties
        self.hairpinPenalty = 0

        # Hyperparameters
        if hyperparams:
            self.hairpinHyperparams = hyperparams.hairpin

    ### Homopolymer Penalty ###

    ### Hairpin Penalty ###

    def hairpin_penalty(self, curElem, elems1, elems2, isKey = False):
        hairpinStats = self.forward_hairpin_counter(curElem, elems1, elems2, isKey) + self.backward_hairpin_counter(curElem, elems1, elems2, isKey)
        return hairpinStats

    def is_in_curElem(self, info, pos):
        return pos >= 0 and pos < info['elem1Size']

    def is_in_elem2(self, info, pos): # is in key
        return (pos % self.motifSize) >= info['elem1Size']

    def is_in_other_elem1(self, info, pos):
        return not (self.is_in_curElem(info, pos) or self.is_in_elem2(info, pos))

    def is_in_same_elem1(self, info, pos1, pos2):
        # there can be 3 possible elem1 in the same check window
        if self.is_in_elem2(info, pos1) or self.is_in_elem2(info, pos2):
            return False
        is_in_same_curElem = self.is_in_curElem(info, pos1) and self.is_in_curElem(info, pos2)
        pos1_is_in_back_elem1 = (not self.is_in_curElem(info, pos1)) and pos1 < 0
        pos2_is_in_back_elem1 = pos2 < 0 and pos2 >= - info['elem2Size'] - info['elem1Size']
        pos1_is_in_front_elem1 = (not self.is_in_curElem(info, pos1)) and pos1 >= info['elem1Size']
        pos2_is_in_front_elem1 = pos2 > info['elem1Size'] and pos1 < info['elem2Size'] + info['elem1Size'] * 2
        return is_in_same_curElem or (pos1_is_in_back_elem1 and pos2_is_in_back_elem1) \
                or (pos1_is_in_front_elem1 and pos2_is_in_front_elem1)

    def is_in_same_elem2(self, info, pos1, pos2):
        # there can only be 2 elem2 in the same check window
        if not (self.is_in_elem2(info, pos1) and self.is_in_elem2(info, pos2)):
            return False
        pos1_is_in_back_elem2 = pos1 < 0
        pos2_is_in_back_elem2 = pos2 < 0 and pos2 >= - info['elem2Size']
        pos1_is_in_front_elem2 = pos1 >= info['elem1Size']
        pos2_is_in_front_elem2 = pos2 > 0 and pos2 < info['elem1Size'] + info['elem2Size']
        return (pos1_is_in_back_elem2 and pos2_is_in_back_elem2) \
                or (pos1_is_in_front_elem2 and pos2_is_in_front_elem2)
    
    # to be used if stuck or start
    def send_to_all_check(self, info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE1_1="", curE1_2="", curE2_1="", curE2_2=""):
        if curJ < 0 and hairpinLength > 0:
            if info['hairpinCount']:
                hairpins.append(1)
            else:
                hairpins.append(self.hairpinHyperparams.hyperparameter**(hairpinLength / self.maxHairpin))

        for j in range(curJ, -1, -1):

            if (not self.send_to_all_pos1_elem1(info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE1_1, curE1_2, curE2_2)) \
                and (not self.send_to_all_pos1_elem2(info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE2_1, curE2_2, curE1_2)):
                curJ -= 1
            else:
                return True

            if curJ < 0 and hairpinLength > 0:
                if info['hairpinCount']:
                    hairpins.append(1)
                else:
                    hairpins.append(self.hairpinHyperparams.hyperparameter**(hairpinLength / self.maxHairpin))
                return True

        return False

    def send_to_all_pos1_elem1(self, info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE1_1="", curE1_2="", curE2=""):
        return self.send_to_elem1_elem1(info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE1_1, curE1_2) or \
            self.send_to_elem1_elem2(info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE1_1, curE2)

    def send_to_elem1_elem1(self, info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE1_1="", curE1_2=""):
        stem1Pos = stem1Start + curJ
        stem2Pos = stem2Start + self.maxHairpin - 1 - curJ

        if self.is_in_elem2(info, stem1Pos) or self.is_in_elem2(info, stem2Pos):
            return False

        if curE1_1 and curE1_2:
            self.hairpin_count_elem1_elem1(info, curE1_1, curE1_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        elif curE1_1:
            if not self.is_in_curElem(info, stem2Pos):
                for e1 in info['elems1']:
                    if self.is_in_same_elem1(info, stem1Pos, stem2Pos) and curE1_1 != e1:
                        continue
                    self.hairpin_count_elem1_elem1(info, curE1_1, e1, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
            self.hairpin_count_elem1_elem1(info, curE1_1, info['curElem'], curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        elif curE1_2:
            if not self.is_in_curElem(info, stem1Pos):
                for e1 in info['elems1']:
                    if self.is_in_same_elem1(info, stem1Pos, stem2Pos) and e1 != curE1_2:
                        continue
                    self.hairpin_count_elem1_elem1(info, e1, curE1_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
            self.hairpin_count_elem1_elem1(info, info['curElem'], curE1_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        else:
            p1_p2_in_curElem = self.is_in_curElem(info, stem1Pos) and self.is_in_curElem(info, stem2Pos)
            if (not p1_p2_in_curElem) and self.is_in_curElem(info, stem1Pos):
                for e1_2 in info['elems1']:
                    if self.is_in_same_elem1(info, stem1Pos, stem2Pos) and info['curElem'] != e1_2:
                        continue
                    self.hairpin_count_elem1_elem1(info, info['curElem'], e1_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
            elif (not p1_p2_in_curElem) and self.is_in_curElem(info, stem2Pos):
                for e1_1 in info['elems1']:
                    if self.is_in_same_elem1(info, stem1Pos, stem2Pos) and info['curElem'] != e1_1:
                        continue
                    self.hairpin_count_elem1_elem1(info, e1_1, info['curElem'], curJ, stem1Start, stem2Start, hairpinLength, hairpins)
            elif not (self.is_in_curElem(info, stem1Pos) or self.is_in_curElem(info, stem2Pos)):
                for e1_1 in info['elems1']:
                    for e1_2 in info['elems1']:
                        if self.is_in_same_elem1(info, stem1Pos, stem2Pos) and e1_1 != e1_2:
                            continue
                        self.hairpin_count_elem1_elem1(info, e1_1, e1_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
            self.hairpin_count_elem1_elem1(info, info['curElem'], info['curElem'], curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        return True

    def send_to_elem1_elem2(self, info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE1="", curE2=""):
        stem1Pos = stem1Start + curJ
        stem2Pos = stem2Start + self.maxHairpin - 1 - curJ

        if not ((not self.is_in_elem2(info, stem1Pos)) and self.is_in_elem2(info, stem2Pos)):
            return False

        if curE1 and curE2:
            self.hairpin_count_elem1_elem2(info, curE1, curE2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        elif curE1:
            for e2 in info['elems2']:
                self.hairpin_count_elem1_elem2(info, curE1, e2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        elif curE2:
            if not self.is_in_curElem(info, stem1Pos):
                for e1 in info['elems1']:
                    self.hairpin_count_elem1_elem2(info, e1, curE2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
            self.hairpin_count_elem1_elem2(info, info['curElem'], curE2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        else:
            for e2 in info['elems2']:
                if not self.is_in_curElem(info, stem1Pos):
                    for e1 in info['elems1']:
                        self.hairpin_count_elem1_elem2(info, e1, e2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
                self.hairpin_count_elem1_elem2(info, info['curElem'], e2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        
        if not (curE2 or info['elems2']):
            return self.send_to_all_check(info, curJ-1, stem1Start, stem2Start, hairpinLength, hairpins, curE1_1=curE1, curE2_2=curE2)
        
        return True

    def send_to_all_pos1_elem2(self, info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE2_1="", curE2_2="", curE1=""):
        return self.send_to_elem2_elem1(info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE2_1, curE1) or \
            self.send_to_elem2_elem2(info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE2_1, curE2_2)

    def send_to_elem2_elem1(self, info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE2="", curE1=""):
        stem1Pos = stem1Start + curJ
        stem2Pos = stem2Start + self.maxHairpin - 1 - curJ
        if not (self.is_in_elem2(info, stem1Pos) and self.is_in_other_elem1(info, stem2Pos)):
            return False
        if curE2 and curE1:
            self.hairpin_count_elem2_elem1(info, curE2, curE1, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        elif curE2:
            if not self.is_in_curElem(info, stem2Pos):
                for e1 in info['elems1']:
                    self.hairpin_count_elem2_elem1(info, curE2, e1, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
            self.hairpin_count_elem2_elem1(info, curE2, info['curElem'], curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        elif curE1:
            for e2 in info['elems2']:
                self.hairpin_count_elem2_elem1(info, e2, curE1, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        else:
            for e2 in info['elems2']:
                if not self.is_in_curElem(info, stem2Pos):
                    for e1 in info['elems1']:
                        self.hairpin_count_elem2_elem1(info, e2, e1, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
                self.hairpin_count_elem2_elem1(info, e2, info['curElem'], curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        
        if not (curE2 or info['elems2']):
            return self.send_to_all_check(info, curJ-1, stem1Start, stem2Start, hairpinLength, hairpins, curE1_2=curE1, curE2_1=curE2)
        
        return True

    def send_to_elem2_elem2(self, info, curJ, stem1Start, stem2Start, hairpinLength, hairpins, curE2_1="", curE2_2=""):
        stem1Pos = stem1Start + curJ
        stem2Pos = stem2Start + self.maxHairpin - 1 - curJ
        if not (self.is_in_elem2(info, stem1Pos) and self.is_in_elem2(info, stem2Pos)):
            return False

        if curE2_1 and curE2_2:
            self.hairpin_count_elem2_elem2(info, curE2_1, curE2_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        elif curE2_1:
            for e2_2 in info['elems2']:
                if self.is_in_same_elem2(info, stem1Pos, stem2Pos) and curE2_1 != e2_2:
                    continue
                self.hairpin_count_elem2_elem2(info, curE2_1, e2_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        elif curE2_2:
            for e2_1 in info['elems2']:
                if self.is_in_same_elem2(info, stem1Pos, stem2Pos) and e2_1 != curE2_2:
                    continue
                self.hairpin_count_elem2_elem2(info, e2_1, curE2_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        else:
            for e2_1 in info['elems2']:
                for e2_2 in info['elems2']:
                    if self.is_in_same_elem2(info, stem1Pos, stem2Pos) and e2_1 != e2_2:
                        continue
                    self.hairpin_count_elem2_elem2(info, e2_1, e2_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins)
        
        if not ((curE2_1 and curE2_2) or info['elems2']):
            return self.send_to_all_check(info, curJ-1, stem1Start, stem2Start, hairpinLength, hairpins, curE2_1=curE2_1, curE2_2=curE2_2)
        
        return True

    def hairpin_count_elem1_elem1(self, info, elem1_1, elem1_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins):
        for j in range(curJ, -1, -1):
            if self.send_to_elem2_elem1(info, j, stem1Start, stem2Start, hairpinLength, hairpins, curE1=elem1_2) or \
            self.send_to_elem2_elem2(info, j, stem1Start, stem2Start, hairpinLength, hairpins) or \
            self.send_to_elem1_elem2(info, j, stem1Start, stem2Start, hairpinLength, hairpins, curE1=elem1_1):
                return

            stem1StartPos = stem1Start % self.motifSize
            stem2StartPos = stem2Start % self.motifSize

            curStem1Pos = (stem1StartPos + j) % self.motifSize
            curStem2Pos = (stem2StartPos + self.maxHairpin - 1 - j) % self.motifSize

            if curStem1Pos >= len(elem1_1) or curStem2Pos >= len(elem1_2):
                continue
            
            if elem1_1[curStem1Pos] == converse[elem1_2[curStem2Pos]]:
                hairpinLength += 1
                if j == 0:
                    if info['hairpinCount']:
                        hairpins.append(1)
                    else:
                        hairpins.append(self.hairpinHyperparams.hyperparameter**(hairpinLength / self.maxHairpin))
            else:
                hairpinLength = 0
                break

    def hairpin_count_elem1_elem2(self, info, elem1, elem2, curJ, stem1Start, stem2Start, hairpinLength, hairpins):
        for j in range(curJ, -1, -1):
            if self.send_to_elem2_elem1(info, j, stem1Start, stem2Start, hairpinLength, hairpins) or \
            self.send_to_elem2_elem2(info, j, stem1Start, stem2Start, hairpinLength, hairpins, curE2_2=elem2) or \
            self.send_to_elem1_elem1(info, j, stem1Start, stem2Start, hairpinLength, hairpins, curE1_1=elem1):
                return

            stem1StartPos = stem1Start % self.motifSize
            stem2StartPos = stem2Start % self.motifSize

            curStem1Pos = (stem1StartPos + j) % self.motifSize
            curStem2Pos = (stem2StartPos + self.maxHairpin - 1 - j) % self.motifSize

            if curStem1Pos >= len(elem1):
                continue
            
            if elem1[curStem1Pos] == converse[elem2[curStem2Pos - info['elem1Size']]]:
                hairpinLength += 1
                if j == 0:
                    if info['hairpinCount']:
                        hairpins.append(1)
                    else:
                        hairpins.append(self.hairpinHyperparams.hyperparameter**(hairpinLength / self.maxHairpin))
            else:
                hairpinLength = 0
                break

    def hairpin_count_elem2_elem1(self, info, elem2, elem1, curJ, stem1Start, stem2Start, hairpinLength, hairpins):
        for j in range(curJ, -1, -1):
            if self.send_to_elem1_elem1(info, j, stem1Start, stem2Start, hairpinLength, hairpins, curE1_2=elem1) or \
            self.send_to_elem1_elem2(info, j, stem1Start, stem2Start, hairpinLength, hairpins) or \
            self.send_to_elem2_elem2(info, j, stem1Start, stem2Start, hairpinLength, hairpins, curE2_1=elem2):
                return

            stem1StartPos = stem1Start % self.motifSize
            stem2StartPos = stem2Start % self.motifSize

            curStem1Pos = (stem1StartPos + j) % self.motifSize
            curStem2Pos = (stem2StartPos + self.maxHairpin - 1 - j) % self.motifSize

            if curStem2Pos >= len(elem1):
                continue
            
            if elem2[curStem1Pos - info['elem1Size']] == converse[elem1[curStem2Pos]]:
                hairpinLength += 1
                if j == 0:
                    if info['hairpinCount']:
                        hairpins.append(1)
                    else:
                        hairpins.append(self.hairpinHyperparams.hyperparameter**(hairpinLength / self.maxHairpin))
            else:
                hairpinLength = 0
                break

    def hairpin_count_elem2_elem2(self, info, elem2_1, elem2_2, curJ, stem1Start, stem2Start, hairpinLength, hairpins):
        for j in range(curJ, -1, -1):
            if self.send_to_elem1_elem1(info, j, stem1Start, stem2Start, hairpinLength, hairpins) or \
            self.send_to_elem1_elem2(info, j, stem1Start, stem2Start, hairpinLength, hairpins, curE2=elem2_2) or \
            self.send_to_elem2_elem1(info, j, stem1Start, stem2Start, hairpinLength, hairpins, curE2=elem2_1):
                return

            stem1StartPos = stem1Start % self.motifSize
            stem2StartPos = stem2Start % self.motifSize

            curStem1Pos = (stem1StartPos + j) % self.motifSize
            curStem2Pos = (stem2StartPos + self.maxHairpin - 1 - j) % self.motifSize
            
            if elem2_1[curStem1Pos - info['elem1Size']] == converse[elem2_2[curStem2Pos - info['elem1Size']]]:
                hairpinLength += 1
                if j == 0:
                    if info['hairpinCount']:
                        hairpins.append(1)
                    else:
                        hairpins.append(self.hairpinHyperparams.hyperparameter**(hairpinLength / self.maxHairpin))
            else:
                hairpinLength = 0
                break

    def backward_hairpin_counter(self, curElem, elems1, elems2, isKey=False, hairpinCount=False):
        hairpins = []
        for i in range(self.maxHairpin):
            stem1Start = len(curElem) - 1 - i
            self.backward_hairpin_counter_at_startPos(curElem, elems1, elems2, stem1Start, hairpins, isKey=isKey, hairpinCount=hairpinCount)

        return np.sum(np.array(hairpins))

    def backward_hairpin_counter_at_startPos(self, curElem, elems1, elems2, stem1Start, hairpins, isKey=False, hairpinCount=True, loopSizeMin=-1, loopSizeMax=-1):
        if loopSizeMin < 0 or loopSizeMax < 0:
            loopSizeMin = self.loopSizeMin
            loopSizeMax = self.loopSizeMax

        for loopSize in range(loopSizeMin, loopSizeMax + 1):
            stem2Start = stem1Start - loopSize - self.maxHairpin
            
            # elem1 has same type (key/payload) as curElem
            # elem2 has opposite type (key/payload) of curElem
            elem1Size = self.keySize if isKey else self.payloadSize
            elem2Size = self.payloadSize if isKey else self.keySize

            info = {'curElem': curElem, 
                    'elems1': elems1, 
                    'elem1Size': elem1Size, 
                    'elems2': elems2, 
                    'elem2Size': elem2Size,
                    'isKey': isKey,
                    'hairpinCount': hairpinCount
                    }
            self.send_to_all_check(info, self.maxHairpin - 1, stem1Start, stem2Start, 0, hairpins)

        return np.sum(np.array(hairpins))

    def forward_hairpin_counter_at_startPos(self, curElem, elems1, elems2, stem1Start, hairpins, isKey=False, hairpinCount=True, loopSizeMin=-1, loopSizeMax=-1):
        if loopSizeMin < 0 or loopSizeMax < 0:
            loopSizeMin = self.loopSizeMin
            loopSizeMax = self.loopSizeMax

        for loopSize in range(loopSizeMin, loopSizeMax + 1):
            stem2Start = stem1Start + loopSize + self.maxHairpin
            
            # elem1 has same type (key/payload) as curElem
            # elem2 has opposite type (key/payload) of curElem
            elem1Size = self.keySize if isKey else self.payloadSize
            elem2Size = self.payloadSize if isKey else self.keySize

            info = {'curElem': curElem, 
                    'elems1': elems1, 
                    'elem1Size': elem1Size, 
                    'elems2': elems2, 
                    'elem2Size': elem2Size,
                    'isKey': isKey,
                    'hairpinCount': hairpinCount
                    }

            self.send_to_all_check(info, self.maxHairpin - 1, stem1Start, stem2Start, 0, hairpins)

        return np.sum(np.array(hairpins))


    def forward_hairpin_counter(self, curElem, elems1, elems2, isKey=False, hairpinCount=False):
        hairpins = []
        for i in range(self.maxHairpin):
            stem1Start = len(curElem) - 1 - i
            self.forward_hairpin_counter_at_startPos(curElem, elems1, elems2, stem1Start, hairpins, isKey=isKey, hairpinCount=hairpinCount)

        return np.sum(np.array(hairpins))


if __name__ == '__main__':
    basePenalties = BasePenalties()
