import keyMotifBuilder
import asyncio

async def generateMotifs(constraints, payloadNum, payloadSize, keyNum, keySize, maxHairpin, gcContentMinPercentage, gcContentMaxPercentage, maxHomopolymer, loopMin, loopMax):
    
    blob = await asyncio.gather(keyMotifBuilder.buildAnswer(constraints, payloadNum, payloadSize, keyNum, keySize, maxHairpin, gcContentMinPercentage, gcContentMaxPercentage, maxHomopolymer, loopMin, loopMax))
    stringBlob = str(blob)
    if stringBlob[len(stringBlob)-4] == 's':
        setThing = stringBlob[2:len(stringBlob)-9]
        isValid  = False
    else:
        setThing = stringBlob[2:len(stringBlob)-8]
        isValid  = True

    keys, payloads, motifs = setThing.split('}, {')
    keys = keys[2:len(keys) - 1]
    payloads = payloads[1:len(payloads) - 1]
    motifs = motifs[1:len(motifs) - 2]

    finalIsValid = isValid

    finalKeys = set()
    for k in keys.split("', '"):
        finalKeys.add(k)

    finalPayloads = set()
    for p in payloads.split("', '"):
        finalPayloads.add(p)

    finalMotifs = set()
    for m in motifs.split("', '"):
        finalMotifs.add(m)

    return finalKeys, finalPayloads, finalMotifs, finalIsValid
