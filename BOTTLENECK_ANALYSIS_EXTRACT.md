# Performance Bottleneck Analysis: chira_extract.py

## Executive Summary

This document identifies performance bottlenecks in `chira_extract.py` and suggests optimization strategies. The script processes chimeric read alignments and extracts chimeras/singletons, which can be slow for large input files.

## Critical Bottlenecks

### 1. O(n²) Complexity in Alignment Pair Generation (Lines 212-312)
**Location**: `extract_and_write()` function, line 212
**Issue**: `itertools.combinations(l_read_alignments, 2)` generates all pairs of alignments for each read
- **Complexity**: O(n²) where n = number of alignments per read
- **Impact**: For reads with many alignments (e.g., 20 alignments = 190 pairs), this becomes expensive
- **Current**: Creates all pairs upfront, then filters
- **Optimization**: Early filtering before pair generation, or use generator instead of list

**Code**:
```python
alignment_pairs = list(itertools.combinations(l_read_alignments, 2))  # O(n²)
for alignment1, alignment2 in alignment_pairs:
    # ... processing ...
```

**Suggested Optimization**:
- Use generator instead of list: `itertools.combinations(l_read_alignments, 2)` (no `list()`)
- Pre-filter alignments before generating pairs (remove invalid ones early)
- Early exit conditions (e.g., skip if prob == 0 before pair generation)

---

### 2. Repeated String Operations (Lines 214-250, 317-319)
**Location**: `extract_and_write()` function
**Issue**: Multiple `.split()`, `.rstrip()` calls on the same strings
- **Impact**: String operations are expensive, especially when repeated
- **Current**: Each alignment line is split multiple times
- **Optimization**: Cache split results, use pre-parsed data structures

**Code**:
```python
# Line 216: Split once
[segmentid1, transcriptid1, ...] = alignment1.rstrip('\n').split('\t')
# Line 244-250: Split again if switch_alignments
[segmentid1, transcriptid1, ...] = alignment2.rstrip('\n').split('\t')
```

**Suggested Optimization**:
- Parse each alignment line once into a dictionary or named tuple
- Cache parsed alignments to avoid re-parsing
- Pre-strip and split all alignments at the start of the function

---

### 3. Linear Search in `guess_region()` Function (Lines 41-118)
**Location**: `guess_region()` function
**Issue**: Nested loops searching through UTR/CDS/mature miRNA annotations
- **Complexity**: O(m) where m = number of UTR/CDS features per transcript
- **Impact**: Called for every alignment, can be thousands of times
- **Current**: Linear search through lists of features
- **Optimization**: Use interval trees or sorted lists with binary search

**Code**:
```python
if transcriptid in d_transcript_annotations['UTR']:
    for pos in d_transcript_annotations['UTR'][transcriptid]:  # Linear search
        # ... check overlap ...
```

**Suggested Optimization**:
- Build interval trees for UTR/CDS features per transcript (using `intervaltree` package)
- Or sort features by start position and use binary search
- Cache region lookups in `d_regions` dictionary (already partially done)

---

### 4. Repeated `extract_annotations()` Calls (Lines 257-264, 321-324)
**Location**: `extract_and_write()` function
**Issue**: `extract_annotations()` called multiple times per read, with repeated dictionary lookups
- **Impact**: Dictionary lookups and `guess_region()` calls are expensive
- **Current**: Called separately for each alignment pair
- **Optimization**: Cache annotation results, batch lookups

**Code**:
```python
geneid1, name1, region1, tx_len1 = extract_annotations(transcriptid1, genomic_pos1, d_regions, f_gtf)
geneid2, name2, region2, tx_len2 = extract_annotations(transcriptid2, genomic_pos2, d_regions, f_gtf)
```

**Suggested Optimization**:
- Pre-extract annotations for all unique (transcriptid, genomic_pos) pairs in the read
- Cache results in a dictionary keyed by `(transcriptid, genomic_pos)`
- Reuse cached results for repeated transcript/position combinations

---

### 5. FASTA Parsing with SeqIO (Line 421)
**Location**: `hybridize_and_write()` function
**Issue**: Uses `SeqIO.parse()` which is slower than raw file parsing
- **Impact**: Biopython's SeqIO is convenient but slower than raw parsing
- **Current**: `SeqIO.parse()` for reading FASTA files
- **Optimization**: Replace with raw file parsing (like in `chira_collapse.py`)

**Code**:
```python
fa_seq = SeqIO.parse(open(os.path.join(outdir, "loci.fa.") + n, 'r', buffering=BUFFER_SIZE), 'fasta')
for record in fa_seq:
    # ... process ...
```

**Suggested Optimization**:
- Parse FASTA manually: read lines, detect headers (starting with `>`), extract sequences
- Similar to optimization in `chira_collapse.py` (2-5x faster)
- Use f-strings for string formatting

---

### 6. No Batch Writing (Lines 348, 352)
**Location**: `extract_and_write()` function
**Issue**: Writing line by line instead of batching
- **Impact**: System calls for each write operation
- **Current**: `fh_chimeras.write("\t".join(a) + "\n")` for each chimera
- **Optimization**: Collect lines in buffer, write in batches

**Code**:
```python
for a in l_best_chimeras:
    fh_chimeras.write("\t".join(a) + "\n")  # One write per line
```

**Suggested Optimization**:
- Collect output lines in a list (e.g., `output_buffer`)
- Write in batches (e.g., every 10,000 lines) using `writelines()`
- Similar to optimization in `chira_quantify.py`

---

### 7. Repeated Float Conversions (Lines 251-255, 326-327)
**Location**: `extract_and_write()` function
**Issue**: Converting strings to floats multiple times, then formatting
- **Impact**: Float conversion and string formatting are expensive
- **Current**: `float(prob1)`, then `"{:.2f}".format(float(prob1))`
- **Optimization**: Convert once, cache results

**Code**:
```python
first_locus_score = float("{:.2f}".format(float(prob1)))  # Convert twice
```

**Suggested Optimization**:
- Parse alignment fields once into typed values (float, int, str)
- Cache float values to avoid repeated conversions
- Format only when needed for output

---

### 8. File I/O for Large Input Files (Lines 374-398)
**Location**: `write_chimeras()` function
**Issue**: Reading entire CRL file line by line, even for chunk processing
- **Impact**: For very large files, I/O becomes a bottleneck
- **Current**: Reads entire file, skips lines outside chunk range
- **Optimization**: Use memory-mapped files for very large files, or seek to chunk start

**Code**:
```python
for line in fh_crl:
    # Process line, but skip if outside chunk range
    if chunk_start > read_count + 1:
        read_count += 1
        continue
```

**Suggested Optimization**:
- Use memory-mapped files for files > 100MB (similar to `chira_quantify.py`)
- Or pre-index file by readid to enable seeking (more complex)

---

### 9. String Formatting in Output (Lines 296-310, 332-339)
**Location**: `extract_and_write()` function
**Issue**: Multiple string conversions and joins
- **Impact**: String operations add up for large datasets
- **Current**: `str()`, `",".join()`, `"\t".join()` operations
- **Optimization**: Pre-format strings, use f-strings, cache formatted values

**Code**:
```python
chimera = [readid, transcriptid1, transcriptid2,
           ..., str(tx_pos_start1), str(tx_pos_end1), ...,
           ",".join([str(arm1_start), str(arm1_end), ...])]
```

**Suggested Optimization**:
- Pre-format numeric values as strings once
- Use f-strings for complex formatting: `f"{arm1_start},{arm1_end},..."`
- Cache formatted strings if reused

---

### 10. Dictionary Lookups in `extract_annotations()` (Lines 167-188)
**Location**: `extract_annotations()` function
**Issue**: Multiple dictionary lookups with fallbacks
- **Impact**: Dictionary lookups are fast, but repeated lookups add up
- **Current**: Multiple `.get()` calls with fallbacks
- **Optimization**: Use `.get()` with default values, cache lookups

**Code**:
```python
geneid = d_transcript_annotations['gid'][transcriptid] if transcriptid in d_transcript_annotations['gid'] else 'NA'
name = d_gene_annotations['name'][geneid] if geneid in d_gene_annotations['name'] else 'NA'
```

**Suggested Optimization**:
- Use `.get()` method: `d_transcript_annotations['gid'].get(transcriptid, 'NA')`
- More Pythonic and slightly faster
- Cache gene-level annotations if same gene appears multiple times

---

## Optimization Priority

### High Priority (Biggest Impact)
1. **Cache string parsing results** (Bottleneck #2): Parse alignment lines once, reuse
2. **Batch writing** (Bottleneck #6): Write in batches of 10,000 lines
3. **Pre-extract annotations** (Bottleneck #4): Cache annotation results per read
4. **Replace SeqIO with raw parsing** (Bottleneck #5): 2-5x faster FASTA reading

### Medium Priority (Moderate Impact)
5. **Optimize `guess_region()`** (Bottleneck #3): Use interval trees or sorted lists
6. **Cache float conversions** (Bottleneck #7): Parse once, reuse
7. **Use generators for combinations** (Bottleneck #1): Avoid creating full list
8. **String formatting optimizations** (Bottleneck #9): Use f-strings, pre-format

### Low Priority (Smaller Impact, but Easy Wins)
9. **Dictionary lookup optimizations** (Bottleneck #10): Use `.get()` method
10. **Memory-mapped files** (Bottleneck #8): For very large input files

---

## Expected Performance Improvements

### Conservative Estimates (Implementing High Priority Only)
- **2-3x faster** overall processing for typical datasets
- **3-5x faster** FASTA parsing (replacing SeqIO)
- **20-30% faster** I/O operations (batch writing)
- **10-20% faster** string operations (caching parsing)

### Aggressive Estimates (Implementing All Optimizations)
- **5-10x faster** overall processing for large datasets
- **10-20x faster** for reads with many alignments (O(n²) optimization)
- **50-70% reduction** in memory usage (generators, caching)

---

## Implementation Notes

### 1. Alignment Parsing Cache
```python
def parse_alignment_line(line):
    """Parse alignment line once, return dict."""
    fields = line.rstrip('\n').split('\t')
    return {
        'segmentid': fields[0],
        'transcriptid': fields[1],
        'locusid': fields[2],
        'crlid': fields[3],
        'tx_pos_start': fields[4],
        'tx_pos_end': fields[5],
        'tx_pos_strand': fields[6],
        'cigar': fields[7],
        'genomic_pos': fields[8],
        'locuspos': fields[9],
        'locusshare': float(fields[10]),
        'prob': float(fields[11]),
        'tpm': float(fields[12])
    }
```

### 2. Batch Writing
```python
BATCH_SIZE = 10000
output_buffer = []
for a in l_best_chimeras:
    output_buffer.append("\t".join(a) + "\n")
    if len(output_buffer) >= BATCH_SIZE:
        fh_chimeras.writelines(output_buffer)
        output_buffer = []
if output_buffer:
    fh_chimeras.writelines(output_buffer)
```

### 3. Raw FASTA Parsing
```python
# Replace SeqIO.parse with raw parsing
with open(fasta_file, 'r', buffering=BUFFER_SIZE) as fh:
    current_id = None
    current_seq = []
    for line in fh:
        if line.startswith('>'):
            if current_id:
                d_loci_seqs[current_id] = ''.join(current_seq).upper().replace('T', 'U')
            current_id = line[1:].strip()[:-3] if len(line) > 4 else line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_id:
        d_loci_seqs[current_id] = ''.join(current_seq).upper().replace('T', 'U')
```

### 4. Annotation Caching
```python
# Pre-extract annotations for all unique (transcriptid, genomic_pos) in read
annotation_cache = {}
for alignment in l_read_alignments:
    parsed = parse_alignment_line(alignment)
    key = (parsed['transcriptid'], parsed['genomic_pos'])
    if key not in annotation_cache:
        annotation_cache[key] = extract_annotations(
            parsed['transcriptid'], parsed['genomic_pos'], d_regions, f_gtf
        )
```

---

## Testing Recommendations

1. **Profile the code** using `cProfile` to identify actual bottlenecks
2. **Test with realistic data sizes**: Small (1K reads), Medium (100K reads), Large (1M+ reads)
3. **Measure improvements** for each optimization separately
4. **Verify correctness**: Ensure optimizations don't change output
5. **Memory profiling**: Check memory usage doesn't increase significantly

---

## References

- Similar optimizations in `chira_quantify.py`: Memory-mapped files, batch writing, string optimizations
- Similar optimizations in `chira_collapse.py`: Raw FASTA parsing instead of SeqIO
- Interval tree library: `pip install intervaltree` (optional, for `guess_region()` optimization)

