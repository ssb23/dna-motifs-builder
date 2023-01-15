import pytest
import basePayloadPenalties as bp
import constraints as c
import hyperparameters as h


def get_constraints(maxHairpin=2, loopSize=1, payloadSize=10, keySize=2, payloadNum=5, maxHom=1, minGC=25, maxGC=60):
    payloadSize = payloadSize
    payloadNum = payloadNum
    maxHom = maxHom
    maxHairpin = maxHairpin
    loopSize = loopSize
    minGc = minGC
    maxGc = maxGC
    keySize = keySize
    
    return c.Constraints(payloadSize, payloadNum, maxHom, maxHairpin, loopSize, minGc, maxGc, keySize)

def get_hyperparameters(hairpinHyperparam=5, homHyperparam=5, gcHyperparam=5):
    hyperparams = {'hom': homHyperparam, 'motifGcContent': gcHyperparam, 'hairpin': hairpinHyperparam}
    return h.Hyperparameters(hyperparams)

def get_score(elemLength, elemHyperparam, maxElem):
    if elemLength == 0:
        return 0
    return elemHyperparam**(elemLength/maxElem)

###### GC content tests ######

def get_score_gc_content(gcCount, curMotifSize, minGc, maxGc, motifSize):
    weight = curMotifSize / motifSize
    minGcContent = (100 * gcCount) / curMotifSize
    maxGcContent = (100 * gcCount) / curMotifSize
    score = 0
    score = max(score, weight * (minGc - minGcContent))
    score = max(score, weight * (maxGcContent - maxGc))
    return score

@pytest.mark.asyncio
async def test_gc_within_bounds():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'AAGG', 'AACC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GA'
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = 0
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_empty_returns_keys_gc_min():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'AAAA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = ''
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = get_score_gc_content(0, 4, minGC, maxGC, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_empty_returns_keys_gc_max():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'AAGC', 'GGCC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = ''
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = get_score_gc_content(4, 4, minGC, maxGC, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_empty_returns_keys_gc_within_bounds():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'AGTG', 'ACCT'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = ''
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = 0
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_edge_min():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 10
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'ACGA', 'TGTC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'AAAAAA'
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = get_score_gc_content(2, 10, minGC, maxGC, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_edge_max():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 10
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'ACGA', 'TGTC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'AGGAGC'
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = get_score_gc_content(6, 10, minGC, maxGC, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_bounds_full_payload():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 6
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'ACGA', 'TGTC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'AGGAGA'
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = get_score_gc_content(5, 10, minGC, maxGC, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_min_full_payload():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 6
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'AAGA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'AAAAAA'
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = get_score_gc_content(1, 10, minGC, maxGC, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_max_full_payload():
    minGC = 20
    maxGC = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 6
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'GGAG', 'AAAA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GGGGGC'
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = get_score_gc_content(9, 10, minGC, maxGC, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_min_and_max_full_payload():
    minGC = 40
    maxGC = 50
    gcHyperparam = 5
    keySize = 4
    payloadSize = 6
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGC=minGC, maxGC=maxGC, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    joints = {'GGGG', 'AAAA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'AAAAGG'
    gcScore = basePayloadPenalties.motif_GC_content_penalty(curPayload)
    result = get_score_gc_content(2, 10, minGC, maxGC, motifSize)
    assert result == gcScore

###### Homopolymer tests ######

@pytest.mark.asyncio
async def test_hom_size_one():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AATT', 'GGCC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(1, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_isolates_long():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AATT', 'GGCC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GAAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(4, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_with_end_joint():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AAGG', 'GGCC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GAAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(6, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_with_full_key_start():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AAAA', 'GGCC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GAAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(8, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_with_full_key_start_with_itself():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AAAA', 'GGCC'}
    payloads = {'TGAAAA', 'CGAAAA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'AGAAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(9, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_with_full_key_start_with_itself_and_other_payloads():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AAAA', 'GGCC'}
    payloads = {'AGAATA', 'AGAACA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'AGAAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(9, homHyperparam, maxHom) * 3
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_full_motif():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AAAA', 'GGCC'}
    payloads = {'AAAATA', 'TGAACA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'AAAAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(10 * maxHom, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_with_full_key():
    maxHom = 2
    homHyperparam = 5
    keySize = 6
    payloadSize = 6
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AAAAAA', 'GGCCCC'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'AAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(10, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_with_full_key_and_with_other_payloads():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AAAA', 'GGCC'}
    payloads = {'TGAATA', 'TGACAA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'AAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(8, homHyperparam, maxHom) + get_score(9, homHyperparam, maxHom) + get_score(10, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_with_full_payload_and_part_key():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'TAAA', 'GGCC'}
    payloads = {'TGAATA', 'TGACAA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'AAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(7, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_with_full_key_and_with_full_other_payloads():
    maxHom = 2
    homHyperparam = 5
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    joints = {'AAAA', 'GGCC'}
    payloads = {'AAAAAA', 'TGACTT'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'AAAA'
    homScore = basePayloadPenalties.homopolymer_penalty(curPayload)
    result = get_score(14, homHyperparam, maxHom) + get_score(8, homHyperparam, maxHom)
    assert result == homScore

##### Hairpin tests ######

##### Backward hairpin tests ######

@pytest.mark.asyncio
async def test_stem1_in_payload_stem2_in_key_edge():
    maxHairpin = 2
    loopSize = 1
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GA'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload_stem2_non_existent():
    maxHairpin = 2
    loopSize = 5
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GA'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(0, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_empty_curPayload_returns_zero():
    maxHairpin = 2
    loopSize = 1
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = ''
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(0, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_empty_key_hairpin_between_payloads():
    maxHairpin = 2
    loopSize = 5
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = set()
    payloads = {'AAAAAAAAAA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'TTT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(1, hairpinHyperparam, maxHairpin) + get_score(2, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload_stem2_in_payload_whole():
    maxHairpin = 4
    loopSize = 2
    payloadSize = 4
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'TT', 'AC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GATC'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(4, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_full_hairpin_stem1_in_payload_stem2_in_key():
    maxHairpin = 2
    loopSize = 1
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GAT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload_stem2_in_payload_and_key_with_2_keys_danger():
    maxHairpin = 2
    loopSize = 1
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'AAT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(1, hairpinHyperparam, maxHairpin) * 2 + get_score(2, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload_stem2_in_key_and_payload_part_hairpin_with_2_keys_danger():
    maxHairpin = 2
    loopSize = 1
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'ATT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(1, hairpinHyperparam, maxHairpin) * 2
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload1_stem2_in_payload1():
    maxHairpin = 2
    loopSize = 1
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'ATTAT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin) 
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_full_hairpin_in_payload_at_edge():
    maxHairpin = 2
    loopSize = 1
    payloadSize = 5
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'ATTAT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin) 
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload1_stem2_in_payload2_part_in_key():
    maxHairpin = 2
    loopSize = 2
    payloadSize = 5
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = {'ATCGG', 'ATCCG'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'CT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(1, hairpinHyperparam, maxHairpin) * (len(payloads) + 1)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload1_and_key_stem2_in_payload2():
    maxHairpin = 3
    loopSize = 2
    payloadSize = 8
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = {'TATAAGGA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'CT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(3, hairpinHyperparam, maxHairpin) + get_score(1, hairpinHyperparam, maxHairpin) * len(joints)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload1_stem2_in_key():
    maxHairpin = 2
    loopSize = 1
    payloadSize = 8
    hairpinHyperparam = 5
    keySize = 4
    constraints = get_constraints(maxHairpin, loopSize, payloadSize, keySize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AAGT', 'ATGC'}
    payloads = {'TATAAGGA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'CT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload1_and_key_stem2_in_key():
    maxHairpin = 3
    loopSize = 1
    payloadSize = 8
    hairpinHyperparam = 5
    keySize = 6
    constraints = get_constraints(maxHairpin, loopSize, payloadSize, keySize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'TAGATT', 'TAGGGC'}
    payloads = {'TATAAGGA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'CT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(3, hairpinHyperparam, maxHairpin) * 2
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload1_stem2_in_payload1_and_payload2(): # impossible
    maxHairpin = 4
    loopSize = 1
    payloadSize = 7
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'CC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'ATAATAT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(4, hairpinHyperparam, maxHairpin) * 2
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_payload1_stem2_in_payload1_and_payload2():
    maxHairpin = 4
    loopSize = 1
    payloadSize = 7
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AA', 'GC'}
    payloads = {'TCATCGA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'TGATTT'
    hairpinScore = basePayloadPenalties.backward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(4, hairpinHyperparam, maxHairpin) + get_score(3, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

###### Forward hairpin tests ######

@pytest.mark.asyncio
async def test_forward_empty_curPayload_returns_zero():
    maxHairpin = 2
    loopSize = 1
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = ''
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(0, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_stem1_in_payload_stem2_in_payload_whole():
    maxHairpin = 4
    loopSize = 2
    payloadSize = 4
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'TT', 'AC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GATC'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(4, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_full_hairpin_stem1_in_payload_stem2_in_key():
    maxHairpin = 2
    loopSize = 7
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'GAT'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_stem1_in_payload_stem2_in_payload_and_key():
    maxHairpin = 2
    loopSize = 8
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AA', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'TAT'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin) + get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_stem1_in_payload_stem2_in_payload_and_key_with_2_keys_danger():
    maxHairpin = 2
    loopSize = 6 + 12 * 2
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'TG', 'GC'}
    payloads = {'AATCATCGAA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'TAT'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin) + get_score(1, hairpinHyperparam, maxHairpin) 
    assert result == hairpinScore


@pytest.mark.asyncio
async def test_forward_stem1_in_payload1_stem2_in_payload1():
    maxHairpin = 2
    loopSize = 7
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'ATTAT'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin) 
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_full_hairpin_in_payload_at_edge():
    maxHairpin = 2
    loopSize = 2
    payloadSize = 5
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'ATTAT'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin) 
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_stem1_in_payload1_and_key_stem2_in_payload2():
    maxHairpin = 3
    loopSize = 2
    payloadSize = 5
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'GC'}
    payloads = {'TTAGA'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'TTTCT'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(3, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_stem1_in_payload1_stem2_in_key():
    maxHairpin = 2
    loopSize = 1
    payloadSize = 8
    hairpinHyperparam = 5
    keySize = 4
    constraints = get_constraints(maxHairpin, loopSize, payloadSize, keySize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AAGT', 'ATAC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'TATAAGTA'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(2, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_stem1_in_payload1_and_key_stem2_in_key():
    maxHairpin = 3
    loopSize = 1
    payloadSize = 8
    hairpinHyperparam = 5
    keySize = 6
    constraints = get_constraints(maxHairpin, loopSize, payloadSize, keySize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'TAATCT', 'TAGGGC'}
    payloads = set()
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    curPayload = 'TATAAGGA'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(3, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_stem1_in_payload1_stem2_in_payload12_and_payload2():
    maxHairpin = 4
    loopSize = 8
    payloadSize = 7
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, payloadSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    joints = {'AT', 'CC'}
    payloads = {'CCATCGG'}
    basePayloadPenalties = bp.BasePayloadPenalties(constraints, hyperparams)
    await basePayloadPenalties.add_joints(joints)
    await basePayloadPenalties.add_payloads(payloads)
    curPayload = 'CTAGATC'
    hairpinScore = basePayloadPenalties.forward_hairpin_counter(curPayload, payloads, joints)
    result = get_score(4, hairpinHyperparam, maxHairpin) * 2
    assert result == hairpinScore
