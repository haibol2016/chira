#!/usr/bin/env python
import argparse
import re
import subprocess
from Bio import SeqIO
from datetime import datetime

# Compile regex patterns once for better performance
_CIGAR_PATTERN_QUERY = re.compile(r'(\d+)([HIMSX=])')
_CIGAR_PATTERN_ALIGN = re.compile(r'(\d+)([DIMX=])')
_CIGAR_PATTERN_END = re.compile(r'(\d+)([DMNX=])')

# Cache for get_bedtools_command to avoid repeated subprocess calls
_bedtools_command_cache = {}


def score_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def overlap(f, s):
    return max(0, min(f[1], s[1]) - max(f[0], s[0]))


def median(x):
    n = len(x)
    mid = int(n/2)
    if not n % 2:
        return (x[mid-1] + x[mid]) / 2.0
    return x[mid]


def query_length(cigar, is_reverse):
    cigar_tup = _CIGAR_PATTERN_QUERY.findall(cigar)
    read_length = 0
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        read_length += int(c[0])
    return read_length


def match_positions(cigar, is_reverse):
    cigar_tup = _CIGAR_PATTERN_QUERY.findall(cigar)
    match_start = match_end = 0
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        if c[1] == "S" or c[1] == "H":
            if not match_start:
                match_start = int(c[0]) + 1
                match_end = int(c[0])
        elif c[1] == "M":
            if match_start:
                match_end = match_end + int(c[0])
            else:
                match_start = 1
                match_end = int(c[0])
        elif c[1] == "I":
            match_end = match_end + int(c[0])
    return match_start, match_end


def is_chimeric(cigar1, cigar2, is_reverse1, is_reverse2, max_allowed_overlap):
    match_start1, match_end1 = match_positions(cigar1, is_reverse1)
    match_start2, match_end2 = match_positions(cigar2, is_reverse2)
    chimeric = True
    if overlap([match_start1, match_end1], [match_start2, match_end2]) > max_allowed_overlap:
        chimeric = False
    return chimeric


def alignment_length(cigar):
    # everything except clipped
    cigar_tup = _CIGAR_PATTERN_ALIGN.findall(cigar)
    align_len = 0
    for c in cigar_tup:
        align_len += int(c[0])
    return align_len


def alignment_end(start, cigar, is_reverse):
    cigar_tup = _CIGAR_PATTERN_END.findall(cigar)
    end = int(start)
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        end += int(c[0])
    return end - 1


def bedentry(referenceid, reference_start, reference_end, readid, strand, cigarstring):
    line = '\t'.join([referenceid,
                      reference_start,
                      reference_end,
                      ','.join([readid, referenceid, reference_start, reference_end, strand, cigarstring]),
                      "1",
                      strand])
    return line


def reverse_complement(seq):
    tab = str.maketrans("ACTGactg", "TGACtgac")
    return seq.translate(tab)[::-1]


def extract_reflengths(ref_fasta, d_reflen):
    # Use context manager for proper file handling
    with open(ref_fasta, 'r') as fh:
        fa_ref = SeqIO.parse(fh, 'fasta')
        for record in fa_ref:
            d_reflen[record.id] = len(record)
    return


def print_w_time(message):
    # Use f-string for slightly better performance than concatenation
    print(f"[{datetime.now().strftime('%d-%m-%Y %H:%M:%S')}] {message}")


def get_bedtools_command(tool_name):
    """
    Get the appropriate BEDTools command format for different versions.
    Older versions use individual commands (e.g., intersectBed, fastaFromBed),
    newer versions use unified bedtools command (e.g., bedtools intersect, bedtools getfasta).
    
    Args:
        tool_name: Name of the tool (e.g., 'intersect', 'getfasta')
    
    Returns:
        List of command parts to use with subprocess or os.system
    """
    # Check cache first to avoid repeated subprocess calls
    if tool_name in _bedtools_command_cache:
        return _bedtools_command_cache[tool_name]
    
    # Map tool names to old and new command formats
    tool_map = {
        'intersect': ('intersectBed', ['bedtools', 'intersect']),
        'getfasta': ('fastaFromBed', ['bedtools', 'getfasta'])
    }
    
    if tool_name not in tool_map:
        raise ValueError("Unknown BEDTools tool: " + tool_name)
    
    old_cmd, new_cmd = tool_map[tool_name]
    
    # Try new format first (bedtools <tool>)
    try:
        subprocess.run(new_cmd + ['--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=2, check=False)
        result = new_cmd
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        # Fall back to old format (individual command)
        result = [old_cmd]
    
    # Cache the result
    _bedtools_command_cache[tool_name] = result
    return result