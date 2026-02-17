# Bottleneck Analysis: chira_merge.py

## Executive Summary

This document identifies performance bottlenecks in `chira_merge.py` and suggests optimization strategies. The script processes BED files to merge alignments and convert coordinates, with several computationally expensive operations.

## Critical Bottlenecks

### 1. **`filter_alignments()` Function - O(n²) Pair Generation** ⚠️ HIGH PRIORITY

**Location:** Lines 36-137

**Problem:**
- Uses `itertools.product()` and `itertools.combinations()` which create full lists in memory
- For reads with many alignments, this creates O(n²) pairs
- Example: 20 alignments → 190 pairs (combinations) or 400 pairs (product)
- Repeated string splitting: `refpos.split(',')[0]` called multiple times per iteration

**Current Code:**
```python
alignments_pairs = list(itertools.combinations(prev_read_alignments, 2))  # Creates full list
# ...
for refpos1 in prev_read_alignments[readpos1]:
    for refpos2 in prev_read_alignments[readpos2]:
        if refpos1.split(',')[0] not in refids1:  # Repeated splitting
```

**Impact:**
- Memory: O(n²) pairs stored in memory
- CPU: Repeated string splitting operations
- For 1000 reads with 10 alignments each: ~45,000 pairs per read = 45M pairs total

**Optimization Strategy:**
1. **Use generator instead of list**: `itertools.combinations()` directly (lazy evaluation)
2. **Pre-parse refpos**: Extract refid once, cache in dictionary
3. **Early exit conditions**: Skip pairs that can't match before expensive operations

**Expected Speedup:** 2-5x for reads with many alignments

---

### 2. **`write_segments()` Function - Repeated Overlap Calculations** ⚠️ HIGH PRIORITY

**Location:** Lines 140-202

**Problem:**
- Nested loops with sorted operations on every iteration
- Overlap calculations performed multiple times for same pairs
- String operations in inner loops

**Current Code:**
```python
for (current_match_start, current_match_end, current_strand) in sorted(filtered_alignments):
    for segmentid in d_read_segments:
        l_read_pos = sorted(d_read_segments[segmentid])  # Sorted every iteration!
        # ... overlap calculations ...
        overlap = max(0, min(first_match_end, current_match_end) - max(first_match_start, current_match_start))
        # ... more overlap calculations ...
```

**Impact:**
- CPU: O(n log n) sorting on every iteration
- CPU: Repeated overlap calculations
- For 1000 segments: ~1M overlap calculations

**Optimization Strategy:**
1. **Pre-sort segments once**: Sort `d_read_segments[segmentid]` once, not every iteration
2. **Cache overlap results**: Store overlap calculations in dictionary
3. **Early exit**: Skip segments that can't overlap based on position

**Expected Speedup:** 3-10x for reads with many segments

---

### 3. **`filter_alignments()` - Repeated String Splitting** ⚠️ MEDIUM PRIORITY

**Location:** Lines 52, 56, 85, 107, 127

**Problem:**
- `refpos.split(',')[0]` called multiple times for same refpos
- String splitting is expensive (creates new list, allocates memory)

**Current Code:**
```python
if refpos.split(',')[0] in refids1:  # Split #1
    # ...
elif refpos.split(',')[0] in refids2:  # Split #2 (same refpos!)
```

**Impact:**
- CPU: Redundant string operations
- For 1M alignments: ~2-3M redundant splits

**Optimization Strategy:**
1. **Pre-parse refpos**: Extract refid once, store in dictionary
2. **Cache split results**: `refid = refpos.split(',')[0]` once per refpos

**Expected Speedup:** 1.5-2x for alignment filtering

---

### 4. **`reads_to_segments()` - Line-by-Line Processing** ⚠️ MEDIUM PRIORITY

**Location:** Lines 271-300

**Problem:**
- Processes one line at a time
- Calls `filter_alignments()` for each read (can be expensive)
- String operations in loop

**Current Code:**
```python
for line in fh_in:
    desc = line.split("\t")[3].split(",")  # Split twice
    # ... process line ...
    prev_read_alignments[match_start, match_end, strand][referencepos] = line
```

**Impact:**
- I/O: Already optimized with buffer
- CPU: String splitting operations

**Optimization Strategy:**
1. **Cache split results**: Store parsed fields in structured format
2. **Batch processing**: Process multiple reads at once (if possible)

**Expected Speedup:** 1.2-1.5x

---

### 5. **`transcript_to_genomic_pos()` - String Operations** ⚠️ MEDIUM PRIORITY

**Location:** Lines 856-950

**Problem:**
- Line-by-line processing with string operations
- Set operations for junction detection
- String concatenation in loops

**Current Code:**
```python
for line in fh_overlapout:
    f = line.rstrip('\n').split('\t')  # Split every line
    # ... string operations ...
    b = "\t".join([id2chr[transcriptid], str(genomicstart), str(genomicend), readid, "1", pol]) + "\n"
```

**Impact:**
- CPU: String operations for every line
- For 10M lines: ~10M splits and joins

**Optimization Strategy:**
1. **Use f-strings**: Faster than `join()` for small lists
2. **Cache parsed fields**: Store parsed data in structured format
3. **Batch writing**: Collect lines, write in batches

**Expected Speedup:** 1.3-2x

---

### 6. **`merge_overlapping_intervals()` - Sorted Operations** ⚠️ LOW PRIORITY

**Location:** Lines 303-326

**Problem:**
- Sorts positions for every transcript
- Overlap calculations in loop

**Current Code:**
```python
for currentpos in sorted(d_desc[transcript], key=lambda tup: tup[0]):  # Sorted every call
    # ... overlap calculations ...
```

**Impact:**
- CPU: O(n log n) sorting per transcript
- For 10K transcripts: ~10K sorts

**Optimization Strategy:**
1. **Pre-sort transcripts**: Sort once before processing
2. **Cache sorted positions**: Store sorted positions in dictionary

**Expected Speedup:** 1.2-1.5x

---

## Recommended Optimization Priority

### Priority 1 (High Impact, Easy Implementation):
1. ✅ **Pre-parse refpos in `filter_alignments()`**: Extract refid once, cache in dictionary
2. ✅ **Use generator for combinations**: `itertools.combinations()` directly (no `list()`)
3. ✅ **Pre-sort segments in `write_segments()`**: Sort once, not every iteration
4. ✅ **Cache string split results**: Store parsed fields in structured format

### Priority 2 (Medium Impact, Medium Implementation):
5. ✅ **Early exit conditions**: Skip pairs/segments that can't match
6. ✅ **Batch writing**: Collect lines, write in batches (10K lines)
7. ✅ **Use f-strings**: Replace `join()` with f-strings for small lists

### Priority 3 (Lower Impact, Easy Implementation):
8. ✅ **Pre-sort transcripts**: Sort once before processing
9. ✅ **Cache overlap calculations**: Store results in dictionary

---

## Implementation Notes

### 1. Pre-parse refpos
```python
# Before:
if refpos.split(',')[0] in refids1:

# After:
d_refpos_to_refid = {refpos: refpos.split(',')[0] for refpos in all_refpos}
if d_refpos_to_refid[refpos] in refids1:
```

### 2. Use generator for combinations
```python
# Before:
alignments_pairs = list(itertools.combinations(prev_read_alignments, 2))

# After:
alignments_pairs = itertools.combinations(prev_read_alignments, 2)  # Generator
```

### 3. Pre-sort segments
```python
# Before:
for segmentid in d_read_segments:
    l_read_pos = sorted(d_read_segments[segmentid])  # Every iteration

# After:
d_read_segments_sorted = {segid: sorted(positions) for segid, positions in d_read_segments.items()}
for segmentid in d_read_segments_sorted:
    l_read_pos = d_read_segments_sorted[segmentid]  # Pre-sorted
```

### 4. Cache string splits
```python
# Before:
desc = line.split("\t")[3].split(",")

# After:
line_fields = line.rstrip('\n').split('\t')
desc = line_fields[3].split(",")  # Cache line_fields for reuse
```

---

## Expected Overall Performance Improvement

- **Best case (reads with many alignments)**: 5-10x speedup
- **Average case**: 2-4x speedup
- **Worst case (reads with few alignments)**: 1.5-2x speedup

---

## Memory Considerations

- **Current**: O(n²) pairs stored in memory for `filter_alignments()`
- **After optimization**: O(n) memory (generator-based)
- **Memory reduction**: 50-90% for reads with many alignments

---

## Testing Recommendations

1. **Test with reads having many alignments** (10+ alignments per read)
2. **Test with reads having few alignments** (1-3 alignments per read)
3. **Test with large BED files** (100M+ lines)
4. **Verify output is identical** to original (byte-for-byte comparison)

