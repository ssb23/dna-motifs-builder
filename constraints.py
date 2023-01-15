class Constraints:
    def __init__(self, payloadSize=5, payloadNum=1, maxHom=2, maxHairpin=2, loopSize=-1, minGc=20, maxGc=60, keySize=1, keyNum=1, loopSizeMin=-1, loopSizeMax=-1):
        # Payload information
        self.payloadSize = payloadSize
        self.payloadNum = payloadNum
        
        # Motifs information
        self.motifSize = payloadSize + keySize
        self.maxPossibleMotifNum = payloadNum * keyNum * keyNum

        # Keys information
        self.keySize = keySize
        self.keyNum = keyNum

        # Constraints for keys and motifs
        self.maxHom = maxHom
        self.maxHairpin = maxHairpin
        # If give loop size range:
        self.loopSizeMin = loopSizeMin if loopSize == -1 else loopSize
        self.loopSizeMax = loopSizeMax if loopSize == -1 else loopSize
        self.minGc = minGc
        self.maxGc = maxGc

        # Assertions
        assert(maxHairpin <= self.motifSize)
        assert(keySize % 2 == 0)
        assert(keySize > 1)
        assert(keySize < self.motifSize)
        assert(self.loopSizeMin <= self.loopSizeMax)
        assert(self.loopSizeMin >= 0)
        assert(payloadSize > 0)
        assert(keyNum > 0)
        assert(payloadNum > 0)
        assert(minGc <= maxGc)
        assert(minGc >= 0)
        assert(maxGc <= 100)
        

if __name__ == '__main__':
    constraints = Constraints()
    