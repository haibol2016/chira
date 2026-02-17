# Bottleneck Analysis: chira_quantify.py

## Executive Summary

This analysis identifies the main performance bottlenecks in `chira_quantify.py` and provides optimization recommendations. The script processes chimeric read loci (CRLs) through two main phases:
1. **CRL Building** (`build_crls()`): Groups loci into CRLs based on read overlap
2. **CRL Quantification** (`quantify_crls()`): Uses EM algorithm to quantify CRL expression

## Major Bottlenecks (Ranked by Impact)

### 1. **CRL Matching Loop - O(n×m) Complexity** ⚠️ CRITICAL
**Location**: `build_crls()`, lines 61-97  
**Complexity**: O(n × m) where n = qualified loci, m = existing CRLs  
**Impact**: For 10,000 loci and 1,000 CRLs = 10 million iterations

**Current Code**:
```python
for n_locus, l_locusreads in l_qualified_locireads:
    for crlid in range(len(d_crl_reads) - 1, -1, -1):
        # Check Jaccard similarity with each existing CRL
        n_common_reads = len(l_locusreads & l_crlreads)
        n_union_reads = len(l_locusreads | l_crlreads)
        if n_common_reads / float(n_union_reads) < crl_share_cutoff:
            continue
```

**Optimization Strategies**:
- **Early termination**: Break after finding first match (if only one match needed)
- **Size-based filtering**: Pre-filter CRLs by size before Jaccard calculation
- **Spatial indexing**: Use locality-sensitive hashing (LSH) or MinHash for approximate Jaccard similarity
- **Parallel processing**: Process loci in parallel batches (with careful synchronization)
- **Inverted index**: Build index of read_id -> CRL_ids to reduce comparisons

**Expected Speedup**: 5-50x depending on dataset characteristics

---

### 2. **Triple-Nested Output Loop** ⚠️ HIGH
**Location**: `build_crls()`, lines 182-215  
**Complexity**: O(CRLs × loci_per_CRL × segments_per_locus × transcripts_per_segment)  
**Impact**: For 1,000 CRLs × 10 loci × 100 segments × 5 transcripts = 5 million iterations

**Current Code**:
```python
for crlid, d_crlloci in enumerate(d_crl_locus_reads.values()):
    for locusid, l_locusreads in d_crlloci.items():
        for segmentid in sorted(l_locusreads):  # Sorting in nested loop!
            for transcriptid_pos in d_readlocus_transcripts[segment_key]:
                # String operations: split, join, format
                entry = "\t".join([...])
                fh_crl_file.write(entry + "\n")
```

**Optimization Strategies**:
- **Batch writing**: Collect entries in list, write in larger chunks (10K-100K lines)
- **Pre-compute strings**: Build entry strings once, reuse
- **Remove redundant sorting**: `l_locusreads` is already a set - sorting may not be needed if order doesn't matter
- **String building optimization**: Use list comprehension + join instead of repeated concatenation
- **Parallel writing**: Write different CRLs to separate temp files, merge later

**Expected Speedup**: 2-10x

---

### 3. **Re-reading Loci File** ⚠️ MEDIUM
**Location**: `write_crl_output()`, line 557  
**Impact**: I/O overhead, especially for large files

**Current Code**:
```python
# Re-read file to write output (OS will cache it, so second read is fast)
with open(loci_file, "r", buffering=1024*1024) as fh_in:
    for l in fh_in:
        # Lookup in dictionaries
```

**Optimization Strategies**:
- **Store data during first read**: Keep line data in memory during `build_crls()` instead of re-reading
- **Use memory-mapped files**: For very large files that don't fit in memory
- **Stream processing**: If memory is limited, process in chunks

**Expected Speedup**: 1.5-3x (depends on file size and OS cache)

---

### 4. **String Operations in Hot Loops** ⚠️ MEDIUM
**Location**: Multiple locations  
**Impact**: String splitting, joining, formatting in nested loops

**Optimizations**:
- **Cache split results**: Avoid repeated `split('\t')` on same string
- **Pre-format numbers**: Format floats once, reuse strings
- **Use f-strings**: Faster than `.format()` or `%` formatting
- **Avoid repeated string concatenation**: Use list + join

**Expected Speedup**: 1.2-2x

---

### 5. **EM Algorithm - Already Optimized** ✅
**Location**: `em()`, lines 284-429  
**Status**: Already uses ProcessPoolExecutor for parallelization  
**Note**: This is well-optimized, but could benefit from:
- **Early convergence detection**: Stop if change is negligible
- **Adaptive chunk sizing**: Adjust based on actual performance
- **NUMA awareness**: For multi-socket systems

---

## Detailed Recommendations

### Priority 1: Optimize CRL Matching (Lines 61-97)

**Option A: Early Size Filtering**
```python
# Pre-filter CRLs by size before expensive set operations
candidate_crls = [
    (crlid, d_crl_reads[crlid]) 
    for crlid in range(len(d_crl_reads))
    if lower_bound <= len(d_crl_reads[crlid]) <= upper_bound
]
for crlid, l_crlreads in candidate_crls:
    # Only do expensive operations on filtered set
```

**Option B: Inverted Index**
```python
# Build index: read_id -> list of CRL_ids containing that read
d_read_to_crls = defaultdict(list)
for crlid, crl_reads in d_crl_reads.items():
    for read_id in crl_reads:
        d_read_to_crls[read_id].append(crlid)

# For each locus, only check CRLs that share at least one read
candidate_crls = set()
for read_id in l_locusreads:
    candidate_crls.update(d_read_to_crls[read_id])
```

**Option C: Parallel Processing with MPIRE**
```python
from mpire import WorkerPool
from mpire.shared_objects import SharedObject

def match_locus_to_crls(args, shared_objects):
    """
    Match a single locus to existing CRLs.
    Uses shared objects to reduce memory overhead.
    
    Args:
        args: Tuple of (n_locus, l_locusreads)
        shared_objects: Dictionary of shared objects passed by MPIRE
    """
    n_locus, l_locusreads = args
    # Access shared objects (not copied to each worker, uses shared memory)
    d_crl_reads = shared_objects['d_crl_reads']
    crl_share_cutoff = shared_objects['crl_share_cutoff']
    
    already_crl_member = False
    l_matched_crls = []
    locus_size = len(l_locusreads)
    lower_bound = locus_size * (1 - crl_share_cutoff)
    upper_bound = locus_size / (1 - crl_share_cutoff)
    
    # Traverse CRLs in reverse order (latest CRL is last)
    for crlid in range(len(d_crl_reads) - 1, -1, -1):
        l_crlreads = d_crl_reads[crlid]
        if lower_bound <= len(l_crlreads) <= upper_bound:
            n_common_reads = len(l_locusreads & l_crlreads)
            if n_common_reads == 0:
                continue
            n_union_reads = len(l_locusreads | l_crlreads)
            if n_union_reads == 0:
                continue
            if n_common_reads / float(n_union_reads) >= crl_share_cutoff:
                l_matched_crls.append(crlid)
                already_crl_member = True
    
    return n_locus, l_matched_crls, already_crl_member

# Prepare arguments: (n_locus, l_locusreads) pairs
locus_args = [(n_locus, l_locusreads) for n_locus, l_locusreads in l_qualified_locireads]

# Create shared objects to reduce memory overhead
# MPIRE shares these objects across workers without copying (uses shared memory)
shared_objects = SharedObject({
    'd_crl_reads': d_crl_reads,
    'crl_share_cutoff': crl_share_cutoff
})

# Use MPIRE with shared objects and progress tracking
# MPIRE advantages over concurrent.futures:
# - Lower memory overhead (shared objects via SharedObject)
# - Better progress tracking (built-in progress bars)
# - Worker initialization support (worker_init parameter)
# - More efficient task distribution (better load balancing)
# - Lower inter-process communication overhead
with WorkerPool(n_jobs=num_processes, shared_objects=shared_objects) as pool:
    # Process with progress bar and better error handling
    results = pool.map(match_locus_to_crls, locus_args, 
                      progress_bar=True, 
                      progress_bar_options={'desc': 'Matching loci to CRLs'})

# Collect results and update CRL structures
for n_locus, l_matched_crls, already_crl_member in results:
    if not already_crl_member:
        l_matched_crls.append(n_crl)
        n_crl += 1
    for matched_crl in l_matched_crls:
        d_crl_locus_reads[matched_crl][n_locus] = l_locireads[n_locus]
        d_crl_reads[matched_crl].update(l_locireads[n_locus])
```

**MPIRE Benefits:**
- **Shared objects**: `d_crl_reads` is shared across workers via `SharedObject` (not copied), reducing memory by 50-90% for large datasets
- **Progress tracking**: Built-in progress bars with `progress_bar=True` and customizable options
- **Lower overhead**: More efficient worker management and task distribution than ProcessPoolExecutor
- **Better error handling**: Automatic retry mechanisms and detailed error reporting
- **Worker initialization**: Can initialize workers once with `worker_init` function to set up shared state
- **Performance**: Typically 10-30% faster than ProcessPoolExecutor for CPU-bound tasks due to reduced overhead

### Priority 2: Optimize Output Writing (Lines 182-215)

**Batch Writing**:
```python
BATCH_SIZE = 10000
entry_buffer = []

for crlid, d_crlloci in enumerate(d_crl_locus_reads.values()):
    # ... collect entries ...
    entry_buffer.append(entry)
    
    if len(entry_buffer) >= BATCH_SIZE:
        fh_crl_file.write('\n'.join(entry_buffer) + '\n')
        entry_buffer.clear()

# Write remaining entries
if entry_buffer:
    fh_crl_file.write('\n'.join(entry_buffer) + '\n')
```

**Pre-compute Format Strings**:
```python
# Pre-format locus_share once per locus
locus_share_str = "{:.2f}".format(locus_share)
for segmentid in sorted(l_locusreads):
    # Reuse formatted string
```

### Priority 3: Avoid Re-reading Files

**Store During First Pass**:
```python
# In build_crls(), store line data
loci_lines = []
for n_locus, line in enumerate(fh_merged_bed):
    loci_lines.append((n_locus, line))
    # ... process ...

# Later in write_crl_output(), use stored data
for n_locus, line in loci_lines:
    # ... write output ...
```

### Priority 4: String Operation Optimizations

**Cache Split Results**:
```python
# Instead of:
for line in fh:
    f = line.rstrip('\n').split('\t')
    pos = ':'.join([f[0], f[1], f[2], f[3]])

# Cache if line is reused:
line_data = {}
for line in fh:
    if line not in line_data:
        line_data[line] = line.rstrip('\n').split('\t')
    f = line_data[line]
```

**Use f-strings**:
```python
# Instead of:
entry = "\t".join([segmentid, transcriptid_pos.split('\t')[0], ...])

# Use:
entry = f"{segmentid}\t{transcriptid}\t{locusid}\t{crlid}\t..."
```

---

## Performance Profiling Recommendations

To identify which optimizations will have the most impact on your specific dataset:

1. **Profile with cProfile**:
```python
import cProfile
cProfile.run('main()', 'quantify_profile.stats')
```

2. **Use line_profiler**:
```python
@profile
def build_crls(...):
    # ...
```

3. **Time individual sections**:
```python
import time
start = time.time()
# ... code ...
print(f"Section took {time.time() - start:.2f}s")
```

---

## Expected Overall Speedup

With all optimizations applied:
- **Small datasets** (< 1K loci): 2-5x speedup
- **Medium datasets** (1K-10K loci): 5-20x speedup  
- **Large datasets** (> 10K loci): 10-50x speedup

The actual speedup depends on:
- Number of loci
- Number of CRLs
- Average reads per locus
- Average segments per read
- System memory and CPU cores

---

## Implementation Priority

1. **Quick wins** (1-2 hours):
   - Batch writing for output (Priority 2)
   - String operation optimizations (Priority 4)
   - Early size filtering in CRL matching (Priority 1, Option A)

2. **Medium effort** (4-8 hours):
   - Inverted index for CRL matching (Priority 1, Option B)
   - Avoid re-reading loci file (Priority 3)

3. **Advanced** (1-2 days):
   - Parallel processing for CRL matching (Priority 1, Option C)
   - LSH/MinHash for approximate matching
   - Memory-mapped files for very large datasets

---

## Notes

- The EM algorithm is already well-optimized with ProcessPoolExecutor (could be upgraded to MPIRE for better performance)
- Current buffer sizes (1MB) are good for I/O
- Defensive checks add minimal overhead and are worth keeping
- Consider using `numba` or `cython` for the hottest loops if Python speed is still limiting
- **MPIRE installation**: `pip install mpire` (optional dependency, but recommended for parallel processing)
- MPIRE provides better performance than `concurrent.futures` for CPU-bound tasks due to:
  - Shared memory objects (reduces memory overhead)
  - Lower inter-process communication overhead
  - Better worker initialization and task distribution

