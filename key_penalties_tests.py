import pytest
import baseKeyPenalties as bp
import constraints as c
import hyperparameters as h


def get_constraints(maxHairpin=2, loopSize=1, payloadSize=10, keySize=2, maxHom=1, minGc=25, maxGc=60):
    payloadSize = payloadSize
    payloadNum = 5
    maxHom = maxHom
    maxHairpin = maxHairpin
    loopSize = loopSize
    minGc = minGc
    maxGc = maxGc
    keySize = keySize
    
    return c.Constraints(payloadSize, payloadNum, maxHom, maxHairpin, loopSize, minGc, maxGc, keySize)

def get_hyperparameters(hairpinHyperparam=5, homHyperparam=5, gcHyperparam=5):
    hyperparams = {'hom': homHyperparam, 'keyGcContent': gcHyperparam, 'hairpin': hairpinHyperparam}
    return h.Hyperparameters(hyperparams)

def get_score(elemLength, elemHyperparam, maxElem):
    if elemLength == 0:
        return 0
    return elemHyperparam**(elemLength/maxElem)

###### GC content tests ######

def get_score_gc_content(gcCount, curSize, minGc, maxGc, totalSize):
    weight = curSize / totalSize
    minGcContent = (100 * gcCount) / curSize
    maxGcContent = (100 * gcCount) / curSize
    score = 0
    if minGcContent < minGc:
        score = max(score, weight * (minGc - minGcContent)) # instead of minGC, do with average (center)?
    if maxGcContent > maxGc:
        score = max(score, weight * (maxGcContent - maxGc)) # instead of maxGC, do with average (center)?
    return score

###### GC content within key tests ######

@pytest.mark.asyncio
async def test_gc_within_key_within_bounds():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AGTG', 'GACT'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'G'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_key(curKey + 'A')
    result = get_score_gc_content(1, 2, minGc, maxGc, keySize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_key_min():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AGTG', 'GACT'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'A'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_key(curKey + 'A')
    result = get_score_gc_content(0, 2, minGc, maxGc, keySize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_key_min_edge():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 12
    payloadSize = 5
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'GGGGGGAAAAAA', 'CCCCCCAAAAAA'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'AAACAAAA'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_key(curKey + 'G')
    result = get_score_gc_content(2, 10, minGc, maxGc, keySize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_key_max():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AGTG', 'GACT'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'G'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_key(curKey + 'C')
    result = get_score_gc_content(2, 2, minGc, maxGc, keySize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_key_max_edge():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 12
    payloadSize = 5
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AGAGAGTGTGTG', 'GAGAGACTCTCT'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'CAGCAACA'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_key(curKey + 'G')
    result = get_score_gc_content(6, 10, minGc, maxGc, keySize)
    assert result == gcScore

###### GC content within motif (complementary) tests ######

@pytest.mark.asyncio
async def test_gc_within_motif_within_bounds():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AGTG', 'GACT'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'AG'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_motif(curKey + 'C')
    result = get_score_gc_content(2, 3 + payloadSize, minGc, maxGc, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_motif_min():
    minGc = 30
    maxGc = 80
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AAGG', 'GACT'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'AG'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_motif(curKey + 'A')
    result = get_score_gc_content(0, 3 + payloadSize, minGc, maxGc, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_motif_min_edge():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 12
    payloadSize = 1
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AAAAAAGGGGGG', 'AAAAAACCCCCC'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'ATATATGG'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_motif(curKey + 'A')
    result = get_score_gc_content(2, 10, minGc, maxGc, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_motif_max():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 4
    payloadSize = 5
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AGTG', 'GACC'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = ''
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_motif(curKey + 'G')
    result = get_score_gc_content(3, 3 + payloadSize, minGc, maxGc, motifSize)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_within_motif_max_edge():
    minGc = 20
    maxGc = 60
    gcHyperparam = 5
    keySize = 12
    payloadSize = 1
    motifSize = payloadSize + keySize
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    hyperparams = get_hyperparameters(gcHyperparam=gcHyperparam)
    keys = {'AAAAAAGGAGGG', 'TTTTTTCCATGG'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'AG'
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    gcScore = baseKeyPenalties.key_GC_content_within_motif(curKey + 'T')
    result = get_score_gc_content(6, 10, minGc, maxGc, motifSize)
    assert result == gcScore


###### Homopoylmer tests ######

@pytest.mark.asyncio
async def test_hom_empty_key_returns_zero():
    maxHom = 2
    homHyperparam = 5
    constraints = get_constraints(maxHom=maxHom)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    curKey = ''
    homScore = baseKeyPenalties.homopolymer_penalty(curKey)
    result = get_score(0, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_size_one():
    maxHom = 2
    homHyperparam = 5
    constraints = get_constraints(maxHom=maxHom)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    keys = {'AT', 'GC'}
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'G'
    homScore = baseKeyPenalties.homopolymer_penalty(curKey)
    result = get_score(1, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_full():
    maxHom = 2
    homHyperparam = 5
    keySize = 10
    constraints = get_constraints(maxHom=maxHom, keySize=keySize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    curKey = 'GGGGGGGGG'
    baseKeyPenalties.start_new_key()
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    homScore = baseKeyPenalties.homopolymer_penalty(curKey + 'G')
    result = get_score(10, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_part():
    maxHom = 2
    homHyperparam = 5
    keySize = 10
    constraints = get_constraints(maxHom=maxHom, keySize=keySize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    curKey = 'GGGGTGGGG'
    baseKeyPenalties.start_new_key()
    for b in curKey:
        await baseKeyPenalties.add_base(b)
    homScore = baseKeyPenalties.homopolymer_penalty(curKey + 'G')
    result = get_score(5, homHyperparam, maxHom)
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_empty_returns_zero():
    maxHom = 2
    homHyperparam = 5
    keySize = 10
    constraints = get_constraints(maxHom=maxHom, keySize=keySize)
    hyperparams = get_hyperparameters(homHyperparam=homHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    curKey = ''
    homScore = baseKeyPenalties.homopolymer_penalty(curKey)
    result = 0
    assert result == homScore

###### Hairpin tests ######

###### Backward hairpin tests ######

@pytest.mark.asyncio
async def test_stem1_and_stem2_in_key1_edge():
    maxHairpin = 2
    loopSize = 10
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    keys = {'AT', 'GC'}
    payloads = set()
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'GC'
    hairpinScore = baseKeyPenalties.backward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(2, hairpinHyperparam, maxHairpin) * 2
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_and_stem2_in_key1_part_edge():
    maxHairpin = 2
    loopSize = 11
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    keys = {'AG', 'GC'}
    payloads = set()
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'GC'
    hairpinScore = baseKeyPenalties.backward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_part_stem1_and_stem2_in_key1_part_edge():
    maxHairpin = 2
    loopSize = 10
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    keys = {'AG', 'GC'}
    payloads = set()
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'G'
    hairpinScore = baseKeyPenalties.backward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_empty_stem1_returns_zero():
    maxHairpin = 2
    loopSize = 10
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    keys = {'AG', 'GC'}
    payloads = set()
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = ''
    hairpinScore = baseKeyPenalties.backward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(0, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_stem1_in_key1_and_stem2_in_key2_whole():
    maxHairpin = 2
    loopSize = 16
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize, keySize=6)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    keys = {'ATATAT', 'GCGCGC'}
    payloads = set()
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'GCGCGC'
    hairpinScore = baseKeyPenalties.backward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(2, hairpinHyperparam, maxHairpin) * 2 + get_score(1, hairpinHyperparam, maxHairpin) * 2
    assert result == hairpinScore

###### Forward hairpin tests ######

@pytest.mark.asyncio
async def test_forward_stem1_and_stem2_in_key1_edge():
    maxHairpin = 2
    loopSize = 10
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    keys = {'AT', 'GC'}
    payloads = set()
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'GC'
    hairpinScore = baseKeyPenalties.forward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(2, hairpinHyperparam, maxHairpin) * 2
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_stem1_and_stem2_in_key1_part_edge():
    maxHairpin = 2
    loopSize = 11
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    keys = {'AG', 'GC'}
    payloads = set()
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'GC'
    hairpinScore = baseKeyPenalties.forward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_part_stem1_in_key1_and_stem2_in_key2_part_edge():
    maxHairpin = 2
    loopSize = 10
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    keys = {'AG', 'GC'}
    payloads = set()
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    await baseKeyPenalties.add_keys(keys)
    curKey = 'G'
    hairpinScore = baseKeyPenalties.forward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_part_stem1_and_stem2_in_key1_part_edge():
    maxHairpin = 2
    loopSize = 8
    hairpinHyperparam = 5
    constraints = get_constraints(maxHairpin, loopSize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    keys = set()
    payloads = set()
    curKey = 'GC'
    hairpinScore = baseKeyPenalties.forward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_part_stem1_and_stem2_in_key1_part():
    maxHairpin = 2
    loopSize = 8
    hairpinHyperparam = 5
    keySize = 4
    constraints = get_constraints(maxHairpin, loopSize, keySize=keySize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    keys = set()
    payloads = set()
    curKey = 'GAAC'
    hairpinScore = baseKeyPenalties.forward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_part_stem1_and_full_stem2_in_key1_part():
    maxHairpin = 2
    loopSize = 10
    hairpinHyperparam = 5
    keySize = 6
    constraints = get_constraints(maxHairpin, loopSize, keySize=keySize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    keys = set()
    payloads = set()
    curKey = 'AAGAAC'
    hairpinScore = baseKeyPenalties.forward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_full_stem1_and_full_stem2_in_key1():
    maxHairpin = 2
    loopSize = 12
    hairpinHyperparam = 5
    keySize = 6
    constraints = get_constraints(maxHairpin, loopSize, keySize=keySize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    keys = set()
    payloads = set()
    curKey = 'AAGTAC'
    hairpinScore = baseKeyPenalties.forward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(2, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

@pytest.mark.asyncio
async def test_forward_full_stem1_and_part_stem2_in_key1():
    maxHairpin = 2
    loopSize = 9
    hairpinHyperparam = 5
    keySize = 6
    constraints = get_constraints(maxHairpin, loopSize, keySize=keySize)
    hyperparams = get_hyperparameters(hairpinHyperparam=hairpinHyperparam)
    baseKeyPenalties = bp.BaseKeyPenalties(constraints, hyperparams)
    keys = set()
    payloads = set()
    curKey = 'TAGTAC'
    hairpinScore = baseKeyPenalties.forward_hairpin_counter(curKey, keys, payloads, isKey=True)
    result = get_score(1, hairpinHyperparam, maxHairpin)
    assert result == hairpinScore

