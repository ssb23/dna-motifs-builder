import pytest
import validation as v
import constraints as c

# class TestFunctions(unittest.TestCase):
def get_constraints(maxHairpin=2, loopSize=1, payloadSize=10, keySize=2, payloadNum=5, maxHom=1, minGc=25, maxGc=60):
    payloadSize = payloadSize
    payloadNum = payloadNum
    maxHom = maxHom
    maxHairpin = maxHairpin
    loopSize = loopSize
    minGc = minGc
    maxGc = maxGc
    keySize = keySize
    
    return c.Constraints(payloadSize, payloadNum, maxHom, maxHairpin, loopSize, minGc, maxGc, keySize)

##### Homopolymer tests ######

@pytest.mark.asyncio
async def test_hom_within_bounds():
    maxHom = 5
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    payloads = {'ATATA', 'TATAT'} 
    keys = {'TTAT', 'TATT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    homScore = validate.get_motifs_and_keys_homopolymer_score()
    result = 0
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_within_bounds_edge():
    maxHom = 2
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    payloads = {'TAATA', 'TATAT'} 
    keys = {'TATA', 'TATA'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    homScore = validate.get_motifs_and_keys_homopolymer_score()
    result = 0
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_whole_motif():
    maxHom = 2
    keySize = 4
    payloadSize = 5
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    payloads = {'AAAAA', 'TATAT'} # joints: ATAA
    keys = {'TTAT', 'TTTT'} # motifs: AT AAAAA AA , AA AAAAA AA, ...
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    homScore = validate.get_motifs_and_keys_homopolymer_score()
    result = -constraints.motifSize * constraints.maxPossibleMotifNum
    assert result == homScore

@pytest.mark.asyncio
async def test_hom_whole_motif_and_part_other_motif():
    maxHom = 2
    keySize = 4
    payloadSize = 7
    constraints = get_constraints(maxHom=maxHom, keySize=keySize, payloadSize=payloadSize)
    payloads = {'GCGCGCG', 'ATATATA'}
    keys = {'TATT', 'TTTT'} # motifs: AATA ATATATA AATA, ATATATA AAAA ATATATA, ...
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    homScore = validate.get_motifs_and_keys_homopolymer_score()
    result = -5
    assert result == homScore

##### GC Content tests ######

def get_score_gc_content(gcCount, size, minGc, maxGc):
    minGcContent = (100 * gcCount) / size
    maxGcContent = (100 * gcCount) / size
    score = 0
    score = max(score, minGc - minGcContent)
    score = max(score, maxGcContent - maxGc)
    return -score

### Motifs GC Content tests ###

@pytest.mark.asyncio
async def test_gc_motif_within_bounds():
    minGc = 20
    maxGc = 60
    keySize = 4
    payloadSize = 7
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'GAGAGAG', 'AGAGAGA'}
    keys = {'TATT', 'TTTT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_motifs_gc_score()
    result = 0
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_motif_min():
    minGc = 20
    maxGc = 60
    keySize = 4
    payloadSize = 7
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'AAAAGAG', 'AGAGAGA'}
    keys = {'TATT', 'TTTT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_motifs_gc_score()
    result = get_score_gc_content(2, payloadSize + keySize, minGc, maxGc) * 4
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_motif_min_edge():
    minGc = 20
    maxGc = 60
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'AAAGAG', 'GAGAGA'}
    keys = {'TATT', 'TTTT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_motifs_gc_score()
    result = 0
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_motif_max():
    minGc = 20
    maxGc = 60
    keySize = 4
    payloadSize = 7
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'CCGAGCG', 'AGAGAGA'}
    keys = {'TAGT', 'TTTT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_motifs_gc_score()
    result = get_score_gc_content(7, payloadSize + keySize, minGc, maxGc) * 2
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_motif_max_edge():
    minGc = 20
    maxGc = 60
    keySize = 4
    payloadSize = 6
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'CGAGCG', 'GAGAGA'}
    keys = {'TAGT', 'TTTT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_motifs_gc_score()
    result = 0
    assert result == gcScore

### Key GC Content tests ###

@pytest.mark.asyncio
async def test_gc_key_within_bounds():
    minGc = 20
    maxGc = 60
    keySize = 4
    payloadSize = 7
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'GCGCGCG', 'ATATATA'}
    keys = {'TGTC', 'TCCT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_keys_gc_score()
    result = 0
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_key_min():
    minGc = 20
    maxGc = 60
    keySize = 4
    payloadSize = 7
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'GCGCGCG', 'ATATATA'}
    keys = {'AAAA', 'TCCT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_keys_gc_score()
    result = get_score_gc_content(0, keySize, minGc, maxGc)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_key_min_edge():
    minGc = 20
    maxGc = 60
    keySize = 10
    payloadSize = 7
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'GCGCGCG', 'ATATATA'}
    keys = {'AAAAAAAGAC', 'GAAAAAAAAG'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_keys_gc_score()
    result = 0
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_key_max():
    minGc = 20
    maxGc = 60
    keySize = 4
    payloadSize = 7
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'GCGCGCG', 'ATATATA'}
    keys = {'GCGC', 'TCCT'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_keys_gc_score()
    result = get_score_gc_content(4, keySize, minGc, maxGc)
    assert result == gcScore

@pytest.mark.asyncio
async def test_gc_key_max_edge():
    minGc = 20
    maxGc = 60
    keySize = 10
    payloadSize = 7
    constraints = get_constraints(minGc=minGc, maxGc=maxGc, keySize=keySize, payloadSize=payloadSize)
    payloads = {'GCGCGCG', 'ATATATA'}
    keys = {'CACGAAGGAC', 'GAAAAAAAAG'}
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    gcScore = validate.get_keys_gc_score()
    result = 0
    assert result == gcScore


##### Hairpin tests ######

@pytest.mark.asyncio
async def test_stem1_in_key_stem2_in_key_edge():
    maxHairpin = 2
    loopSize = 10
    hairpinHyperparam = 5
    keySize = 4
    payloadSize = 9
    constraints = get_constraints(maxHairpin, loopSize, keySize=keySize, payloadSize=payloadSize)
    payloads = {'ATAATAAAA', 'AAAAAAATA'}
    keys = {'ATTT', 'TATT', 'TTTT'} # motifs: AATA ATAATAAAA AATA, AAAAAAATA AAAA AAAAAAATA, AAAT AAAAAAATA AAAT, ...
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    hairpinCount = validate.get_motifs_and_keys_hairpin_score()
    result = -3
    assert result == hairpinCount

@pytest.mark.asyncio
async def test_stem1_in_payload_stem2_in_key():
    maxHairpin = 3
    loopSize = 3
    hairpinHyperparam = 5
    keySize = 4
    payloadSize = 8
    constraints = get_constraints(maxHairpin, loopSize, keySize=keySize, payloadSize=payloadSize)
    payloads = {'ATATAGAA', 'AAAAAATA'}
    keys = {'TAGT', 'TATT', 'TTTT'} # motifs: ACTA ATATAGAA AATA, ...
    validate = v.Validate(constraints)
    await validate.add_keys_and_payloads(keys, payloads)
    hairpinCount = validate.get_motifs_and_keys_hairpin_score()
    result = -2
    assert result == hairpinCount
