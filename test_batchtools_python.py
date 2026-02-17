#!/usr/bin/env python
"""
Simple test script for batchtools submission from Python.
Usage: python test_batchtools_python.py
"""
import sys
import os
import tempfile
import json
import subprocess

# Add ChiRA to path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from chira_map import submit_chunks_with_batchtools

class TestArgs:
    """Minimal args object for testing"""
    def __init__(self):
        self.batchtools_queue = "normal"  # CHANGE THIS to your queue name
        self.batchtools_cores = 1
        self.batchtools_memory = "1GB"
        self.batchtools_walltime = "00:05"
        self.batchtools_conda_env = ""
        self.batchtools_max_parallel = 1
        self.batchtools_template = None
        self.outdir = tempfile.mkdtemp(prefix="chira_test_")

def main():
    print("=== Testing Batchtools Submission from Python ===\n")
    
    # Check prerequisites
    print("1. Checking prerequisites...")
    
    # Check R
    try:
        result = subprocess.run(['which', 'Rscript'], capture_output=True, text=True)
        if result.returncode == 0:
            print("   ✓ Rscript found:", result.stdout.strip())
        else:
            print("   ✗ Rscript not found!")
            return 1
    except Exception as e:
        print(f"   ✗ Error checking Rscript: {e}")
        return 1
    
    # Check LSF
    try:
        result = subprocess.run(['which', 'bsub'], capture_output=True, text=True)
        if result.returncode == 0:
            print("   ✓ bsub found:", result.stdout.strip())
        else:
            print("   ✗ bsub not found! LSF may not be configured.")
            return 1
    except Exception as e:
        print(f"   ✗ Error checking bsub: {e}")
        return 1
    
    # Check template
    template_file = os.path.join(script_dir, "lsf_custom.tmpl")
    if os.path.exists(template_file):
        print(f"   ✓ Template file found: {template_file}")
    else:
        print(f"   ⚠ Template file not found: {template_file}")
        print("      Will use built-in 'lsf-simple' template")
    
    print()
    
    # Create test data
    print("2. Creating test data...")
    args = TestArgs()
    
    # Create a dummy chunk file
    chunk_dir = os.path.join(args.outdir, "chunks")
    os.makedirs(chunk_dir, exist_ok=True)
    
    chunk_file = os.path.join(chunk_dir, "chunk_000.fasta")
    with open(chunk_file, 'w') as f:
        f.write(">test_sequence\n")
        f.write("ATCGATCGATCG\n")
    
    chunk_files = [chunk_file]
    alignment_job_types = []  # Empty for testing - we just want to test submission
    per_chunk_processes = 1
    
    print(f"   ✓ Created test chunk: {chunk_file}")
    print(f"   ✓ Output directory: {args.outdir}")
    print()
    
    # Submit jobs
    print("3. Submitting test job via batchtools...")
    print(f"   Queue: {args.batchtools_queue}")
    print(f"   Cores: {args.batchtools_cores}")
    print(f"   Memory: {args.batchtools_memory}")
    print()
    
    try:
        reg_dir, job_ids = submit_chunks_with_batchtools(
            args, chunk_files, chunk_dir, alignment_job_types, per_chunk_processes
        )
        
        if reg_dir and job_ids:
            print(f"   ✓ Submission successful!")
            print(f"   Registry directory: {reg_dir}")
            print(f"   Batchtools job IDs: {job_ids}")
            print()
            
            # Wait a moment for LSF to register
            print("4. Waiting for LSF to register jobs...")
            import time
            time.sleep(3)
            
            # Check LSF
            print("5. Checking LSF for submitted jobs...")
            try:
                result = subprocess.run(
                    ['bjobs', '-u', os.environ.get('USER', '')],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                
                if result.returncode == 0:
                    print("   LSF jobs (bjobs -u $USER):")
                    if result.stdout.strip():
                        print("   " + "\n   ".join(result.stdout.strip().split("\n")))
                    else:
                        print("   (no jobs found)")
                    
                    # Check registry for LSF batch IDs
                    print()
                    print("6. Checking batchtools registry for LSF batch IDs...")
                    job_ids_file = os.path.join(reg_dir, "job_ids.txt")
                    if os.path.exists(job_ids_file):
                        with open(job_ids_file, 'r') as f:
                            bt_job_ids = [line.strip() for line in f if line.strip()]
                        
                        # Try to query batchtools status
                        try:
                            r_script = f'''
library(batchtools)
reg <- loadRegistry("{reg_dir}", writeable=FALSE)
job_table <- getJobTable(reg=reg)
cat(paste(job_table$batch.id, collapse="\\n"))
'''
                            with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
                                f.write(r_script)
                                temp_r = f.name
                            
                            result = subprocess.run(
                                ['Rscript', temp_r],
                                capture_output=True,
                                text=True,
                                timeout=10
                            )
                            os.unlink(temp_r)
                            
                            if result.stdout.strip():
                                lsf_ids = result.stdout.strip().split('\n')
                                print(f"   ✓ Found LSF batch IDs: {', '.join(lsf_ids)}")
                                print(f"   Verify with: bjobs {' '.join(lsf_ids)}")
                            else:
                                print("   ✗ No LSF batch IDs found in registry!")
                                print("   Jobs were not submitted to LSF.")
                        except Exception as e:
                            print(f"   ⚠ Could not query registry: {e}")
                            print(f"   Check registry manually: {reg_dir}")
                else:
                    print(f"   ⚠ bjobs command failed: {result.stderr}")
                    
            except Exception as e:
                print(f"   ⚠ Error checking LSF: {e}")
            
            print()
            print("=== Test Complete ===")
            print(f"Registry directory: {reg_dir}")
            print("Check logs:", os.path.join(reg_dir, "logs"))
            print("Check status in R:")
            print(f'  library(batchtools); reg <- loadRegistry("{reg_dir}"); getStatus(reg)')
            
            return 0
        else:
            print("   ✗ Submission failed!")
            print("   Check error messages above")
            return 1
            
    except Exception as e:
        print(f"   ✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    # Prompt for queue name
    print("NOTE: Make sure to edit the queue name in this script!")
    print("Current queue:", TestArgs().batchtools_queue)
    print()
    
    response = input("Continue with test? (y/n): ")
    if response.lower() != 'y':
        print("Cancelled.")
        sys.exit(0)
    
    sys.exit(main())

