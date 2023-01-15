from motifBuilder import MotifBuilder
from keyBuilder import KeyBuilder
from constraints import Constraints
from validation import Validate
from hyperparameters import Hyperparameters
from language import converse

import numpy as np

def get_constraints(maxHairpin=2, loopSize=1, payloadSize=10, keySize=2, keyNum=5, payloadNum=5, maxHom=1, minGc=25, maxGc=60, loopMin=-1, loopMax=-1):
    payloadSize = payloadSize
    payloadNum = payloadNum
    maxHom = maxHom
    maxHairpin = maxHairpin
    loopMin = loopMin
    loopMax = loopMax
    minGc = minGc
    maxGc = maxGc
    keySize = keySize
    keyNum = keyNum
    
    return Constraints(payloadSize, payloadNum, maxHom, maxHairpin, -1, minGc, maxGc, keySize, keyNum, loopMin, loopMax)

def get_weights(hom=1, keyGcContent=1, motifGcContent=1, hairpin=1):
    weights = {'hom': hom, 'keyGcContent': keyGcContent, 'motifGcContent': motifGcContent, 'hairpin': hairpin}
    return weights

def get_hyperparameters(hom=5, keyGcContent=5, motifGcContent=5, hairpin=5):
    hyperparams = {'hom': hom, 'keyGcContent': keyGcContent, 'motifGcContent': motifGcContent, 'hairpin': hairpin}
    return Hyperparameters(hyperparams)

async def buildKeys(withConstraints, constraints, internalHyperparameters, weights):
    keyBuilder = KeyBuilder(constraints, internalHyperparameters, weights)
    keys = await keyBuilder.buildAllKeys(withConstraints)
    return keys

async def buildMotifs(withConstraints, constraints, internalHyperparameters, weights, joints):
    motifBuilder = MotifBuilder(constraints, internalHyperparameters, weights)
    await motifBuilder.add_joints(joints)
    motifs = await motifBuilder.buildAllPayloads(withConstraints)
    return motifs

def keys_to_joints(keys):
    joints = set()
    for k in keys:
        joint = ''
        for i in range(len(k) - 1, -1, -1):
            joint += converse[k[i]]
        joints.add(joint)
    return joints

async def build_motifs_and_keys(withKeyConstraints, withMotifConstraints, constraints):
    import asyncio

    internalHyperparameters = get_hyperparameters()
    weights = get_weights()
    
    keys = await buildKeys(withKeyConstraints, constraints, internalHyperparameters, weights)
    joints = keys_to_joints(keys)
    payloads = await buildMotifs(withMotifConstraints, constraints, internalHyperparameters, weights, joints)
    
    startJoint = set()
    endJoint = set()
    for j in joints:
        jointSize = int(len(j) / 2)
        endJoint.add(j[:jointSize])
        startJoint.add(j[jointSize:])
    
    motifs = set()
    for sJoint in startJoint:
        for eJoint in endJoint:
            for p in payloads:
                motifs.add(sJoint + p + eJoint)

  #  print('keys: ', keys)
   # print('payloads: ', payloads)
    #print('motifs: ', motifs)
    return keys, payloads, motifs

async def main():
    maxHom = 2
    maxHairpin = 9
    constraints = get_constraints(maxHom=maxHom, maxHairpin=maxHairpin)
    withKeyConstraints = {'hom', 'keyGcContent', 'hairpin'}
    withMotifConstraints = {'hom', 'motifGcContent', 'hairpin'}
    keys, payloads, motifs = await build_motifs_and_keys(withKeyConstraints, withMotifConstraints, constraints)

    # Validation
    keyMotifsValidation = Validate(constraints)
    keyMotifsValidation.add_keys(keys)
    await keyMotifsValidation.add_payloads(payloads)

    homScore = keyMotifsValidation.get_motifs_and_keys_homopolymer_score()
    gcMotifScore = keyMotifsValidation.get_motifs_gc_score()
    gcKeyScore = keyMotifsValidation.get_keys_gc_score()
    hairpinScore = keyMotifsValidation.get_motifs_and_keys_hairpin_score()
    totalScore = homScore + gcKeyScore + gcMotifScore + hairpinScore
    # print('homScore: ', homScore)
    # print('gcMotifScore: ', gcMotifScore)
    # print('hairpinScore: ', hairpinScore)
    # print('gcKeyScore: ', gcKeyScore)
    return  keys, payloads, motifs, totalScore == 0

async def buildAnswer(withConstraints, payloadNum, payloadSize, keyNum, keySize, maxHairpin, gcContentMinPercentage, gcContentMaxPercentage, maxHomopolymer, loopMin, loopMax):
    constraints = get_constraints(maxHairpin=maxHairpin, loopSize=-1, payloadSize=payloadSize, keySize=keySize, keyNum=keyNum, payloadNum=payloadNum, maxHom=maxHomopolymer, minGc=gcContentMinPercentage, maxGc=gcContentMaxPercentage, loopMin=loopMin, loopMax=loopMax)
    withKeyConstraints = set()
    withMotifConstraints = set()
    if 'hom' in withConstraints:
        withKeyConstraints.add('hom')
        withMotifConstraints.add('hom')

    if 'motifGcContent' in withConstraints:
        withKeyConstraints.add('keyGcContent')
        withMotifConstraints.add('motifGcContent')

    if 'hairpin' in withConstraints:
        withKeyConstraints.add('hairpin')
        withMotifConstraints.add('hairpin')

    keys, payloads, motifs = await build_motifs_and_keys(withKeyConstraints, withMotifConstraints, constraints)

    # Validation
    keyMotifsValidation = Validate(constraints)
    keyMotifsValidation.add_keys(keys)
    await keyMotifsValidation.add_payloads(payloads)

    totalScore = 0
    if 'hom' in withConstraints:
        totalScore += keyMotifsValidation.get_motifs_and_keys_homopolymer_score()
    if 'motifGcContent' in withConstraints:
        totalScore += keyMotifsValidation.get_motifs_gc_score()
        totalScore += keyMotifsValidation.get_keys_gc_score()
    if 'hairpin' in withConstraints:
        totalScore += keyMotifsValidation.get_motifs_and_keys_hairpin_score()
    # print('homScore: ', homScore)
    # print('gcMotifScore: ', gcMotifScore)
    # print('hairpinScore: ', hairpinScore)
    # print('gcKeyScore: ', gcKeyScore)
    return  keys, payloads, motifs, totalScore == 0
   

if __name__ == '__main__':
    import asyncio
    asyncio.run(main())

