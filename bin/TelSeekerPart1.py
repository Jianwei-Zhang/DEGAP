#!/usr/bin/env python3
"""
TelSeekerPart1_new.py - Extract telomeric reads using window-based scanning

Modified algorithm:
- Instead of checking fixed left/right ends, scan for telomere markers anywhere in the read
- For each marker found, create a window (500bp upstream + marker + 500bp downstream)
  - If upstream < 500bp, compensate by extending downstream
  - If downstream < 500bp, compensate by extending upstream
  - Target window size: 1000 + 2*len(motif) bp (e.g., 1012bp for TTAGGG)
- Calculate TRC for each window and take the maximum
- Classify based on which marker type (forward/reverse) gave the highest TRC

Input:
    - Reads split into fragments in: --out/hifi_reads_part/ and/or --out/ont_reads_part/
    
Output:
    - --out/part1.telo.reads/left.telo.reads.fa   # Left telomeric reads
    - --out/part1.telo.reads/right.telo.reads.fa  # Right telomeric reads  
    - --out/part1.telo.reads/part1.log            # Processing log

Usage:
    python TelSeekerPart1_new.py --motif TTAGGG --out output_dir -t 100
"""

import os
import sys
import argparse
import glob
import subprocess
import tempfile
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Import TRC calculation functions from TelSeekerCheck
from TelSeekerCheck import calculate_trc


# ===== Motif Filtering Functions =====

def get_reverse_complement(motif: str) -> str:
    """
    Get reverse complement of motif
    
    Args:
        motif: DNA motif
    
    Returns:
        Reverse complement sequence
    """
    return str(Seq(motif).reverse_complement())


def count_exact_motif_matches(sequence: str, motif: str, include_reverse_complement: bool = True) -> int:
    """
    Count exact matches of motif in sequence (non-overlapping)
    Considers both forward motif and its reverse complement
    Returns the maximum count between forward and reverse complement
    
    Args:
        sequence: DNA sequence to search
        motif: Exact motif to search for (e.g., 'TTAGGG')
        include_reverse_complement: Whether to also count reverse complement
    
    Returns:
        Maximum number of non-overlapping exact matches (max of forward and RC)
    """
    sequence = sequence.upper()
    motif = motif.upper()
    
    # Count forward motif (non-overlapping)
    forward_count = sequence.count(motif)
    
    # Count reverse complement if requested and take the maximum
    if include_reverse_complement:
        rc_motif = get_reverse_complement(motif)
        # Only check RC if it's different from forward
        if rc_motif != motif:
            rc_count = sequence.count(rc_motif)
            return max(forward_count, rc_count)
    
    return forward_count


class TeloReadsExtractor:
    """Extract telomeric reads from split read files using window-based scanning"""
    
    def __init__(self, output_dir: str, motif: str, threads: int = 20,
                 trc_threshold: float = 0.7, check_length: int = 1000,
                 batch_size: int = 1000, enable_second_filter: bool = True,
                 min_motif_count: int = 10, overlapping: bool = False):
        """
        Initialize TeloReadsExtractor

        Args:
            output_dir: Output directory (where DEGAP placed reads_part/)
            motif: Telomere repeat motif (e.g., TTAGGG)
            threads: Number of threads for parallel processing (-t parameter)
            trc_threshold: TRC threshold for telomeric classification
            check_length: Length to check at each end (default: 1000bp)
                         Used as window_flank = check_length // 2 for scanning
            batch_size: Number of reads to read at once (default: 1000)
            enable_second_filter: Enable second-level motif count filtering (default: True)
            min_motif_count: Minimum motif count for second filter (default: 10)
            overlapping: Use overlapping count in second filter (default: False)
        """
        self.output_dir = output_dir
        self.motif = motif.upper()
        self.threads = threads
        self.trc_threshold = trc_threshold
        self.check_length = check_length
        # Convert check_length to window_flank (half of check_length)
        self.window_flank = check_length // 2
        self.kmer_length = len(motif) - 2
        self.batch_size = batch_size
        self.enable_second_filter = enable_second_filter
        self.min_motif_count = min_motif_count
        self.overlapping = overlapping
        
        # Generate window markers (2x motif length)
        self.pat_right = self.motif * 2  # e.g., TTAGGGTTAGGG
        self.pat_left = get_reverse_complement(self.pat_right)  # e.g., CCCTAACCCTAA
        
        # Window size = 2 * flank + 2 * motif_length
        self.window_size = 2 * self.window_flank + len(self.pat_right)
        
        # Locate reads_part directories
        self.hifi_reads_part_dir = os.path.join(output_dir, "hifi_reads_part")
        self.ont_reads_part_dir = os.path.join(output_dir, "ont_reads_part")
        
        # Output directory
        self.telo_reads_dir = os.path.join(output_dir, "part1.telo.reads")
        os.makedirs(self.telo_reads_dir, exist_ok=True)
        
        # Output files - Level 1 (temp files)
        self.temp_left_output = os.path.join(self.telo_reads_dir, "temp1.left.telo.reads.fa")
        self.temp_right_output = os.path.join(self.telo_reads_dir, "temp1.right.telo.reads.fa")
        self.temp_left_hifi = os.path.join(self.telo_reads_dir, "temp1.left.hifi.fa")
        self.temp_left_ont = os.path.join(self.telo_reads_dir, "temp1.left.ont.fa")
        self.temp_right_hifi = os.path.join(self.telo_reads_dir, "temp1.right.hifi.fa")
        self.temp_right_ont = os.path.join(self.telo_reads_dir, "temp1.right.ont.fa")
        
        # Output files - Level 2 final outputs
        self.left_output = os.path.join(self.telo_reads_dir, "left.telo.reads.fa")
        self.right_output = os.path.join(self.telo_reads_dir, "right.telo.reads.fa")
        
        # Platform-specific output files (for Global file creation)
        self.left_hifi_output = os.path.join(self.telo_reads_dir, "left.telo.reads.hifi.fa")
        self.left_ont_output = os.path.join(self.telo_reads_dir, "left.telo.reads.ont.fa")
        self.right_hifi_output = os.path.join(self.telo_reads_dir, "right.telo.reads.hifi.fa")
        self.right_ont_output = os.path.join(self.telo_reads_dir, "right.telo.reads.ont.fa")
        
        # Output files - Global merged outputs
        self.global_left_output = os.path.join(self.telo_reads_dir, "Global.left.telo.reads.fa")
        self.global_right_output = os.path.join(self.telo_reads_dir, "Global.right.telo.reads.fa")
        
        # Output files - Level 2 failed (temp files)
        self.temp_left_failed = os.path.join(self.telo_reads_dir, "temp2.left.telo.reads.filter.failed.fa")
        self.temp_right_failed = os.path.join(self.telo_reads_dir, "temp2.right.telo.reads.filter.failed.fa")
        
        self.log_file = os.path.join(self.telo_reads_dir, "part1.log")
        
        # Statistics
        self.stats = {
            'total_reads': 0,
            'left_telo': 0,
            'right_telo': 0,
            'both_telo': 0,
            'non_telo': 0
        }
        
        # Second filter statistics
        self.filter_stats = {
            'left': {'total': 0, 'passed': 0, 'failed': 0},
            'right': {'total': 0, 'passed': 0, 'failed': 0}
        }
        
    def run(self):
        """Run the complete telomeric reads extraction workflow"""
        print(f"[Part1-NEW] Starting telomeric reads extraction...")
        print(f"  Motif: {self.motif}")
        print(f"  Window markers:")
        print(f"    Right (forward): {self.pat_right}")
        print(f"    Left (reverse):  {self.pat_left}")
        print(f"  TRC threshold: {self.trc_threshold}")
        print(f"  Check length: {self.check_length} bp")
        print(f"  Window flank: {self.window_flank} bp (= check_length / 2)")
        print(f"  Window size: {self.window_size} bp (target, with compensation)")
        print(f"  Threads: {self.threads}")
        print(f"  Batch size: {self.batch_size} reads")
        print(f"  Output order: Platform-grouped, reads sorted by ID within each platform")
        print()
        
        # Find all split read files
        split_files = self._find_split_files()
        
        if not split_files:
            print(f"[Part1-NEW] ERROR: No split read files found in:")
            print(f"  {self.hifi_reads_part_dir}")
            print(f"  {self.ont_reads_part_dir}")
            print(f"  Please run DEGAP.py first to create these directories.")
            sys.exit(1)
        
        print(f"[Part1-NEW] Found {len(split_files)} split read files")
        
        # Process split files in parallel
        self._process_split_files_parallel(split_files)
        
        # Merge results
        self._merge_results()
        
        # Apply second-level filter if enabled
        if self.enable_second_filter:
            self._apply_second_filter()
        
        # Create Global merged files
        if self.enable_second_filter:
            self._create_global_files()
        
        # Write log
        self._write_log()
        
        # Print summary
        self._print_summary()
        
        print(f"\n[Part1-NEW] Complete! Results saved to: {self.telo_reads_dir}")
        
    def _find_split_files(self):
        """Find all split read files and categorize by platform"""
        split_files = []
        self.hifi_files = []  # Store HiFi files separately
        self.ont_files = []   # Store ONT files separately
        
        # Check HiFi reads_part
        if os.path.exists(self.hifi_reads_part_dir):
            hifi_files = sorted(glob.glob(os.path.join(self.hifi_reads_part_dir, "*.fa*")))
            self.hifi_files = hifi_files
            split_files.extend(hifi_files)
            print(f"[Part1-NEW] Found {len(hifi_files)} HiFi split files")
        
        # Check ONT reads_part
        if os.path.exists(self.ont_reads_part_dir):
            ont_files = sorted(glob.glob(os.path.join(self.ont_reads_part_dir, "*.fa*")))
            self.ont_files = ont_files
            split_files.extend(ont_files)
            print(f"[Part1-NEW] Found {len(ont_files)} ONT split files")
        
        return split_files
    
    def _process_split_files_parallel(self, split_files):
        """Process split files using GNU parallel"""
        print(f"\n[Part1-NEW] Processing {len(split_files)} files with {self.threads} threads...")

        # Create temporary directory
        temp_dir = os.path.join(self.telo_reads_dir, "temp_parallel")
        os.makedirs(temp_dir, exist_ok=True)

        try:
            # Create processing script
            script_path = self._create_processing_script(temp_dir)

            # Create file list
            file_list_path = os.path.join(temp_dir, "file_list.txt")
            with open(file_list_path, 'w') as f:
                for file_path in split_files:
                    f.write(f"{file_path}\n")

            # Run GNU parallel
            parallel_cmd = [
                'parallel', '--will-cite', '-j', str(self.threads),
                sys.executable, script_path, '{}',
                '::::', file_list_path
            ]
            
            print(f"[Part1-NEW] Running parallel processing...")
            result = subprocess.run(parallel_cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"[Part1-NEW] Warning: parallel processing had issues")
                print(f"  Return code: {result.returncode}")
                if result.stderr:
                    print(f"  Stderr: {result.stderr[:500]}")
            
            print("[Part1-NEW] Parallel processing completed")
            
        finally:
            pass  # Keep temp_dir for merging
    
    def _create_processing_script(self, temp_dir):
        """Create Python script for processing single split file with window-based algorithm"""
        script_path = os.path.join(temp_dir, "process_single_file.py")
        
        # Embed TRC calculation functions directly to avoid import issues
        script_content = f'''#!/usr/bin/env python3
import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq

# ===== Embedded TRC Calculation Functions =====
def _pattern_scramble_telo(pattern, cut_length):
    """Generate all k-mer permutations of telomere motif"""
    pattern = (pattern or "").upper()
    if not pattern:
        return []
    
    cut_length = max(1, cut_length)
    extended = pattern * 2
    kmers = set()
    
    for i in range(len(pattern)):
        fragment = extended[i : i + cut_length]
        if len(fragment) == cut_length:
            kmers.add(fragment)
    
    kmers.update(str(Seq(kmer).reverse_complement()) for kmer in list(kmers))
    return sorted(kmers)

def patterns_to_search_topsicle(motif, cut_length):
    """Generate regex patterns for Topsicle-based telomere detection"""
    kmers = _pattern_scramble_telo(motif, cut_length)
    return [f"(?={{re.escape(kmer)}})" for kmer in kmers]

def calculate_trc(sequence, motif, kmer_length, check_length=None):
    """Calculate Telomere Repeat Content (TRC) using Topsicle method"""
    seq_to_check = sequence.upper()
    if check_length:
        seq_to_check = seq_to_check[:min(len(seq_to_check), check_length)]
    
    if len(seq_to_check) < kmer_length:
        return 0.0
    
    pattern_all = patterns_to_search_topsicle(motif, kmer_length)
    compiled_patterns = [re.compile(pattern) for pattern in pattern_all]
    
    ratio_perfect_hit = max(1.0, len(seq_to_check) / len(motif))
    max_trc = 0.0
    
    for pattern in compiled_patterns:
        matches = len([m.start() for m in pattern.finditer(seq_to_check)])
        pattern_trc = matches / ratio_perfect_hit if ratio_perfect_hit > 0 else 0.0
        max_trc = max(max_trc, pattern_trc)
    
    return max_trc

def find_all_occurrences(seq, pattern):
    """Find all occurrences of pattern in sequence (overlapping allowed for finding positions)"""
    positions = []
    start = 0
    while True:
        pos = seq.find(pattern, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1  # Allow overlapping search
    return positions
# ===== End of Embedded Functions =====

def process_file(input_file, motif, kmer_length, pat_right, pat_left, 
                 window_flank, trc_threshold, temp_dir, batch_size=1000):
    """Process a single split file to extract telomeric reads using window-based scanning"""
    base_name = os.path.basename(input_file)
    left_output = os.path.join(temp_dir, f"left_{{base_name}}")
    right_output = os.path.join(temp_dir, f"right_{{base_name}}")
    stats_output = os.path.join(temp_dir, f"stats_{{base_name}}.txt")
    
    left_reads = []
    right_reads = []
    total = 0
    both_telo = 0
    
    marker_len = len(pat_right)
    
    try:
        # Batch reading: read batch_size reads at a time
        with open(input_file, 'r') as handle:
            record_iter = SeqIO.parse(handle, 'fasta')
            
            while True:
                # Read a batch of reads
                batch = []
                for _ in range(batch_size):
                    try:
                        batch.append(next(record_iter))
                    except StopIteration:
                        break
                
                if not batch:
                    break  # No more reads
                
                # Process each read in the batch
                for record in batch:
                    total += 1
                    seq_str = str(record.seq).upper()
                    
                    # Find all window markers
                    right_positions = find_all_occurrences(seq_str, pat_right)
                    left_positions = find_all_occurrences(seq_str, pat_left)
                    
                    # Calculate TRC for each window
                    max_right_trc = 0.0
                    max_left_trc = 0.0
                    
                    # Target window size: 1000 + 2*len(motif)
                    target_window_size = 2 * window_flank + marker_len
                    
                    # Process right marker windows (forward)
                    for pos in right_positions:
                        # Extract window with fixed size, compensating for insufficient flanks
                        ideal_start = pos - window_flank
                        ideal_end = pos + marker_len + window_flank
                        
                        if ideal_start < 0:
                            # Upstream insufficient, start from 0 and extend downstream
                            window_start = 0
                            window_end = min(len(seq_str), target_window_size)
                        elif ideal_end > len(seq_str):
                            # Downstream insufficient, end at sequence end and extend upstream
                            window_end = len(seq_str)
                            window_start = max(0, len(seq_str) - target_window_size)
                        else:
                            # Both flanks sufficient
                            window_start = ideal_start
                            window_end = ideal_end
                        
                        window_seq = seq_str[window_start:window_end]
                        
                        if len(window_seq) >= marker_len:
                            window_trc = calculate_trc(window_seq, motif, kmer_length)
                            max_right_trc = max(max_right_trc, window_trc)
                    
                    # Process left marker windows (reverse)
                    for pos in left_positions:
                        # Extract window with fixed size, compensating for insufficient flanks
                        ideal_start = pos - window_flank
                        ideal_end = pos + marker_len + window_flank
                        
                        if ideal_start < 0:
                            # Upstream insufficient, start from 0 and extend downstream
                            window_start = 0
                            window_end = min(len(seq_str), target_window_size)
                        elif ideal_end > len(seq_str):
                            # Downstream insufficient, end at sequence end and extend upstream
                            window_end = len(seq_str)
                            window_start = max(0, len(seq_str) - target_window_size)
                        else:
                            # Both flanks sufficient
                            window_start = ideal_start
                            window_end = ideal_end
                        
                        window_seq = seq_str[window_start:window_end]
                        
                        if len(window_seq) >= marker_len:
                            window_trc = calculate_trc(window_seq, motif, kmer_length)
                            max_left_trc = max(max_left_trc, window_trc)
                    
                    # Determine the maximum TRC and its source
                    max_trc = max(max_right_trc, max_left_trc)
                    
                    if max_trc >= trc_threshold:
                        # Determine which type contributed the max TRC
                        if max_right_trc >= trc_threshold and max_left_trc >= trc_threshold:
                            # Both pass threshold, check which is higher
                            if max_right_trc >= max_left_trc:
                                right_reads.append(record)
                            else:
                                left_reads.append(record)
                            # Note: We don't count as "both_telo" anymore since we pick the dominant one
                        elif max_right_trc >= trc_threshold:
                            right_reads.append(record)
                        else:  # max_left_trc >= trc_threshold
                            left_reads.append(record)
        
        # Write results
        with open(left_output, 'w') as f:
            SeqIO.write(left_reads, f, 'fasta')
        
        with open(right_output, 'w') as f:
            SeqIO.write(right_reads, f, 'fasta')
        
        # Write statistics
        with open(stats_output, 'w') as f:
            f.write(f"total={{total}}\\n")
            f.write(f"left={{len(left_reads)}}\\n")
            f.write(f"right={{len(right_reads)}}\\n")
            f.write(f"both={{both_telo}}\\n")
        
        print(f"Processed {{input_file}}: L={{len(left_reads)}}, R={{len(right_reads)}}")
        
    except Exception as e:
        print(f"Error processing {{input_file}}: {{e}}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python process_single_file.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    process_file(
        input_file,
        motif="{self.motif}",
        kmer_length={self.kmer_length},
        pat_right="{self.pat_right}",
        pat_left="{self.pat_left}",
        window_flank={self.window_flank},
        trc_threshold={self.trc_threshold},
        temp_dir="{temp_dir}",
        batch_size={self.batch_size}
    )
'''
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def _merge_results(self):
        """Merge results from all parallel processes, ordered by platform (HiFi first, ONT second)"""
        print(f"\n[Part1-NEW] Merging results...")
        
        temp_dir = os.path.join(self.telo_reads_dir, "temp_parallel")
        
        # Collect statistics
        stats_files = glob.glob(os.path.join(temp_dir, "stats_*.txt"))
        for stats_file in stats_files:
            try:
                with open(stats_file, 'r') as f:
                    for line in f:
                        if '=' in line:
                            key, value = line.strip().split('=')
                            if key == 'total':
                                self.stats['total_reads'] += int(value)
                            elif key == 'left':
                                self.stats['left_telo'] += int(value)
                            elif key == 'right':
                                self.stats['right_telo'] += int(value)
                            elif key == 'both':
                                self.stats['both_telo'] += int(value)
            except Exception as e:
                print(f"Warning: Error reading stats from {stats_file}: {e}")
        
        # Calculate non-telomeric
        self.stats['non_telo'] = (self.stats['total_reads'] - 
                                   self.stats['left_telo'] - 
                                   self.stats['right_telo'] - 
                                   self.stats['both_telo'])
        
        # Merge left telomeric reads to temp1 file (platform-ordered)
        print(f"[Part1-NEW] Merging left telomeric reads (HiFi first, ONT second)...")
        with open(self.temp_left_output, 'w') as outfile, \
             open(self.temp_left_hifi, 'w') as hifi_out, \
             open(self.temp_left_ont, 'w') as ont_out:
            # First: HiFi reads
            for input_file in self.hifi_files:
                base_name = os.path.basename(input_file)
                left_file = os.path.join(temp_dir, f"left_{base_name}")
                if os.path.exists(left_file) and os.path.getsize(left_file) > 0:
                    try:
                        for record in SeqIO.parse(left_file, 'fasta'):
                            SeqIO.write(record, outfile, 'fasta')
                            SeqIO.write(record, hifi_out, 'fasta')  # Also save to HiFi-specific file
                    except Exception as e:
                        print(f"Warning: Error reading {left_file}: {e}")
            
            # Second: ONT reads
            for input_file in self.ont_files:
                base_name = os.path.basename(input_file)
                left_file = os.path.join(temp_dir, f"left_{base_name}")
                if os.path.exists(left_file) and os.path.getsize(left_file) > 0:
                    try:
                        for record in SeqIO.parse(left_file, 'fasta'):
                            SeqIO.write(record, outfile, 'fasta')
                            SeqIO.write(record, ont_out, 'fasta')  # Also save to ONT-specific file
                    except Exception as e:
                        print(f"Warning: Error reading {left_file}: {e}")
        
        # Merge right telomeric reads to temp1 file (platform-ordered)
        print(f"[Part1-NEW] Merging right telomeric reads (HiFi first, ONT second)...")
        with open(self.temp_right_output, 'w') as outfile, \
             open(self.temp_right_hifi, 'w') as hifi_out, \
             open(self.temp_right_ont, 'w') as ont_out:
            # First: HiFi reads
            for input_file in self.hifi_files:
                base_name = os.path.basename(input_file)
                right_file = os.path.join(temp_dir, f"right_{base_name}")
                if os.path.exists(right_file) and os.path.getsize(right_file) > 0:
                    try:
                        for record in SeqIO.parse(right_file, 'fasta'):
                            SeqIO.write(record, outfile, 'fasta')
                            SeqIO.write(record, hifi_out, 'fasta')  # Also save to HiFi-specific file
                    except Exception as e:
                        print(f"Warning: Error reading {right_file}: {e}")
            
            # Second: ONT reads
            for input_file in self.ont_files:
                base_name = os.path.basename(input_file)
                right_file = os.path.join(temp_dir, f"right_{base_name}")
                if os.path.exists(right_file) and os.path.getsize(right_file) > 0:
                    try:
                        for record in SeqIO.parse(right_file, 'fasta'):
                            SeqIO.write(record, outfile, 'fasta')
                            SeqIO.write(record, ont_out, 'fasta')  # Also save to ONT-specific file
                    except Exception as e:
                        print(f"Warning: Error reading {right_file}: {e}")
        
        # Clean up temporary directory
        try:
            shutil.rmtree(temp_dir)
            print(f"[Part1-NEW] Temporary files cleaned up")
        except Exception as e:
            print(f"[Part1-NEW] Warning: Could not remove temp directory: {e}")
    
    def _apply_second_filter(self):
        """Apply second-level motif count filter to telomeric reads"""
        print(f"\n[Part1-NEW] Applying second-level motif count filter...")
        print(f"  Motif: {self.motif}")
        print(f"  Reverse complement: {get_reverse_complement(self.motif)}")
        print(f"  Min motif count: {self.min_motif_count}")
        print(f"  Strategy: max(forward_count, RC_count)")
        print()
        
        # Filter left telomeric reads - check anywhere in sequence
        left_stats = self._filter_reads_by_motif_count(
            input_file=self.temp_left_output,
            output_file=self.left_output,
            failed_file=self.temp_left_failed,
            side='left'
        )
        
        # Filter right telomeric reads - check anywhere in sequence
        right_stats = self._filter_reads_by_motif_count(
            input_file=self.temp_right_output,
            output_file=self.right_output,
            failed_file=self.temp_right_failed,
            side='right'
        )
        
        # Filter platform-specific files (HiFi)
        print(f"\n[Part1-NEW] Filtering platform-specific files...")
        self._filter_reads_by_motif_count(
            input_file=self.temp_left_hifi,
            output_file=self.left_hifi_output,
            failed_file=None,  # Don't save failed for platform files
            side='left-hifi'
        )
        self._filter_reads_by_motif_count(
            input_file=self.temp_right_hifi,
            output_file=self.right_hifi_output,
            failed_file=None,
            side='right-hifi'
        )
        
        # Filter platform-specific files (ONT)
        self._filter_reads_by_motif_count(
            input_file=self.temp_left_ont,
            output_file=self.left_ont_output,
            failed_file=None,
            side='left-ont'
        )
        self._filter_reads_by_motif_count(
            input_file=self.temp_right_ont,
            output_file=self.right_ont_output,
            failed_file=None,
            side='right-ont'
        )
        
        # Store filter statistics
        self.filter_stats['left'] = left_stats
        self.filter_stats['right'] = right_stats
        
        # Clean up temp files
        self._cleanup_temp_files()
        
        print(f"[Part1-NEW] Second-level filtering complete")
    
    def _filter_reads_by_motif_count(self, input_file: str, output_file: str, 
                                      failed_file: str, side: str) -> dict:
        """
        Filter reads based on motif count (anywhere in the full sequence)
        
        Args:
            input_file: Input FASTA file (temp1 files)
            output_file: Output filtered FASTA file (final files)
            failed_file: Output file for failed reads (temp2 files)
            side: 'left' or 'right'
        
        Returns:
            Dictionary with statistics
        """
        if not os.path.exists(input_file):
            return {'total': 0, 'passed': 0, 'failed': 0}
        
        if os.path.getsize(input_file) == 0:
            return {'total': 0, 'passed': 0, 'failed': 0}
        
        print(f"[Filter] Processing {side} telomeric reads...")
        
        total_reads = 0
        passed_reads = 0
        passed_records = []
        failed_records = []
        
        with open(input_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                total_reads += 1
                seq_str = str(record.seq)
                
                # Count motif matches in the ENTIRE sequence (not just ends)
                # This is more consistent with the window-based approach
                motif_count = count_exact_motif_matches(seq_str, self.motif)
                
                # Filter based on threshold
                if motif_count >= self.min_motif_count:
                    passed_reads += 1
                    # Add motif count to description
                    new_record = SeqRecord(
                        seq=record.seq,
                        id=record.id,
                        description=f"{record.description} motif_count={motif_count}"
                    )
                    passed_records.append(new_record)
                else:
                    # Record failed reads
                    failed_record = SeqRecord(
                        seq=record.seq,
                        id=record.id,
                        description=f"{record.description} motif_count={motif_count} FAILED"
                    )
                    failed_records.append(failed_record)
        
        # Sort passed reads by ID for consistent ordering
        passed_records.sort(key=lambda rec: rec.id)
        
        # Write passed reads to final output
        with open(output_file, 'w') as f:
            SeqIO.write(passed_records, f, 'fasta')
        
        # Write failed reads to temp2 file (if specified)
        if failed_file:
            with open(failed_file, 'w') as f:
                SeqIO.write(failed_records, f, 'fasta')
        
        filtered_reads = total_reads - passed_reads
        
        print(f"[Filter] {side.upper()} results:")
        print(f"  Total reads: {total_reads}")
        if total_reads > 0:
            print(f"  Passed (>={self.min_motif_count} motifs): {passed_reads} ({passed_reads/total_reads*100:.1f}%)")
            print(f"  Failed (<{self.min_motif_count} motifs): {filtered_reads} ({filtered_reads/total_reads*100:.1f}%)")
            print(f"  Note: Output sorted by read ID")
        print()
        
        return {
            'total': total_reads,
            'passed': passed_reads,
            'failed': filtered_reads
        }
    
    def _cleanup_temp_files(self):
        """Clean up temporary files after second filter"""
        temp_files = [
            self.temp_left_output,
            self.temp_right_output,
            self.temp_left_failed,
            self.temp_right_failed,
            self.temp_left_hifi,
            self.temp_left_ont,
            self.temp_right_hifi,
            self.temp_right_ont
        ]
        
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except Exception as e:
                    print(f"Warning: Could not remove temp file {temp_file}: {e}")
    
    def _create_global_files(self):
        """
        Create Global merged files with platform grouping:
        - Global.left.telo.reads.fa = HiFi left + RC(HiFi right) + ONT left + RC(ONT right)
        - Global.right.telo.reads.fa = HiFi right + RC(HiFi left) + ONT right + RC(ONT left)
        """
        print(f"\n[Part1-NEW] Creating Global merged files (platform-grouped)...")
        
        # Create Global.left.telo.reads.fa
        print(f"[Global] Creating Global.left.telo.reads.fa...")
        print(f"  Order: HiFi left → RC(HiFi right) → ONT left → RC(ONT right)")
        with open(self.global_left_output, 'w') as out_f:
            # 1. HiFi left reads
            if os.path.exists(self.left_hifi_output):
                for record in SeqIO.parse(self.left_hifi_output, 'fasta'):
                    SeqIO.write(record, out_f, 'fasta')
            
            # 2. RC(HiFi right reads)
            if os.path.exists(self.right_hifi_output):
                for record in SeqIO.parse(self.right_hifi_output, 'fasta'):
                    rc_record = SeqRecord(
                        seq=record.seq.reverse_complement(),
                        id=record.id,
                        description=record.description + " [RC]"
                    )
                    SeqIO.write(rc_record, out_f, 'fasta')
            
            # 3. ONT left reads
            if os.path.exists(self.left_ont_output):
                for record in SeqIO.parse(self.left_ont_output, 'fasta'):
                    SeqIO.write(record, out_f, 'fasta')
            
            # 4. RC(ONT right reads)
            if os.path.exists(self.right_ont_output):
                for record in SeqIO.parse(self.right_ont_output, 'fasta'):
                    rc_record = SeqRecord(
                        seq=record.seq.reverse_complement(),
                        id=record.id,
                        description=record.description + " [RC]"
                    )
                    SeqIO.write(rc_record, out_f, 'fasta')
        
        # Create Global.right.telo.reads.fa
        print(f"[Global] Creating Global.right.telo.reads.fa...")
        print(f"  Order: HiFi right → RC(HiFi left) → ONT right → RC(ONT left)")
        with open(self.global_right_output, 'w') as out_f:
            # 1. HiFi right reads
            if os.path.exists(self.right_hifi_output):
                for record in SeqIO.parse(self.right_hifi_output, 'fasta'):
                    SeqIO.write(record, out_f, 'fasta')
            
            # 2. RC(HiFi left reads)
            if os.path.exists(self.left_hifi_output):
                for record in SeqIO.parse(self.left_hifi_output, 'fasta'):
                    rc_record = SeqRecord(
                        seq=record.seq.reverse_complement(),
                        id=record.id,
                        description=record.description + " [RC]"
                    )
                    SeqIO.write(rc_record, out_f, 'fasta')
            
            # 3. ONT right reads
            if os.path.exists(self.right_ont_output):
                for record in SeqIO.parse(self.right_ont_output, 'fasta'):
                    SeqIO.write(record, out_f, 'fasta')
            
            # 4. RC(ONT left reads)
            if os.path.exists(self.left_ont_output):
                for record in SeqIO.parse(self.left_ont_output, 'fasta'):
                    rc_record = SeqRecord(
                        seq=record.seq.reverse_complement(),
                        id=record.id,
                        description=record.description + " [RC]"
                    )
                    SeqIO.write(rc_record, out_f, 'fasta')
        
        # Count reads in Global files
        global_left_count = sum(1 for _ in SeqIO.parse(self.global_left_output, 'fasta'))
        global_right_count = sum(1 for _ in SeqIO.parse(self.global_right_output, 'fasta'))
        
        print(f"[Global] Global.left.telo.reads.fa: {global_left_count} reads")
        print(f"[Global] Global.right.telo.reads.fa: {global_right_count} reads")
        print(f"[Global] Global files created successfully (platform-grouped)")
    
    def _write_log(self):
        """Write processing log"""
        with open(self.log_file, 'w') as f:
            f.write("TelSeekerPart1_new - Telomeric Reads Extraction Log (Window-based)\n")
            f.write("=" * 60 + "\n\n")
            
            # Algorithm description
            f.write("ALGORITHM: Window-based scanning\n")
            f.write("-" * 60 + "\n")
            f.write(f"Window markers:\n")
            f.write(f"  Right (forward): {self.pat_right}\n")
            f.write(f"  Left (reverse):  {self.pat_left}\n")
            f.write(f"Window configuration:\n")
            f.write(f"  Check length: {self.check_length} bp (total)\n")
            f.write(f"  Flank size: {self.window_flank} bp (= check_length / 2, per side)\n")
            f.write(f"  Target window size: {self.window_size} bp (fixed)\n")
            f.write(f"  Compensation: If one flank is insufficient, extend the other side\n")
            f.write(f"Scan method: Find all markers anywhere in read, calculate TRC for each window\n")
            f.write(f"Classification: Based on maximum TRC across all windows\n")
            f.write(f"Output order:\n")
            f.write(f"  - Platform level: HiFi first, ONT second\n")
            f.write(f"  - Within platform: Sorted by read ID (alphabetical)\n")
            f.write(f"\n")
            
            # First-level TRC filtering parameters and statistics
            f.write("FIRST-LEVEL FILTERING (TRC-based, window scanning)\n")
            f.write("-" * 60 + "\n")
            f.write(f"Parameters:\n")
            f.write(f"  Motif: {self.motif}\n")
            f.write(f"  TRC threshold: {self.trc_threshold}\n")
            f.write(f"  Threads: {self.threads}\n")
            f.write(f"\n")
            f.write(f"Statistics:\n")
            f.write(f"  Total reads processed: {self.stats['total_reads']}\n")
            f.write(f"  Left telomeric reads: {self.stats['left_telo']}\n")
            f.write(f"  Right telomeric reads: {self.stats['right_telo']}\n")
            f.write(f"  Both ends telomeric (excluded): {self.stats['both_telo']}\n")
            f.write(f"  Non-telomeric reads: {self.stats['non_telo']}\n")
            f.write(f"\n")
            f.write(f"Temp output files (Level 1):\n")
            f.write(f"  Left reads: {self.temp_left_output} [cleaned after Level 2]\n")
            f.write(f"  Right reads: {self.temp_right_output} [cleaned after Level 2]\n")
            f.write(f"\n")
            
            # Second-level motif count filtering (if enabled)
            if self.enable_second_filter:
                f.write("SECOND-LEVEL FILTERING (Motif count-based)\n")
                f.write("-" * 60 + "\n")
                f.write(f"Parameters:\n")
                f.write(f"  Motif: {self.motif}\n")
                f.write(f"  Reverse complement: {get_reverse_complement(self.motif)}\n")
                f.write(f"  Min motif count: {self.min_motif_count}\n")
                f.write(f"  Counting method: {'overlapping' if self.overlapping else 'non-overlapping'}\n")
                f.write(f"  Strategy: max(forward_count, RC_count)\n")
                f.write(f"  Scan region: Entire read sequence\n")
                f.write(f"\n")
                
                # Left reads statistics
                left_total = self.filter_stats['left']['total']
                left_passed = self.filter_stats['left']['passed']
                left_failed = self.filter_stats['left']['failed']
                
                f.write(f"Left telomeric reads:\n")
                f.write(f"  Input reads: {left_total}\n")
                if left_total > 0:
                    f.write(f"  Passed (>={self.min_motif_count} motifs): {left_passed} ({left_passed/left_total*100:.2f}%)\n")
                    f.write(f"  Failed (<{self.min_motif_count} motifs): {left_failed} ({left_failed/left_total*100:.2f}%)\n")
                else:
                    f.write(f"  Passed: 0\n")
                    f.write(f"  Failed: 0\n")
                f.write(f"  Input file: {self.temp_left_output}\n")
                f.write(f"  Passed output (final): {self.left_output}\n")
                f.write(f"  Failed output (temp): {self.temp_left_failed} [cleaned]\n")
                f.write(f"\n")
                
                # Right reads statistics
                right_total = self.filter_stats['right']['total']
                right_passed = self.filter_stats['right']['passed']
                right_failed = self.filter_stats['right']['failed']
                
                f.write(f"Right telomeric reads:\n")
                f.write(f"  Input reads: {right_total}\n")
                if right_total > 0:
                    f.write(f"  Passed (>={self.min_motif_count} motifs): {right_passed} ({right_passed/right_total*100:.2f}%)\n")
                    f.write(f"  Failed (<{self.min_motif_count} motifs): {right_failed} ({right_failed/right_total*100:.2f}%)\n")
                else:
                    f.write(f"  Passed: 0\n")
                    f.write(f"  Failed: 0\n")
                f.write(f"  Input file: {self.temp_right_output}\n")
                f.write(f"  Passed output (final): {self.right_output}\n")
                f.write(f"  Failed output (temp): {self.temp_right_failed} [cleaned]\n")
                f.write(f"\n")
                
                # Overall summary
                total_input = left_total + right_total
                total_passed = left_passed + right_passed
                total_failed = left_failed + right_failed
                
                f.write(f"Overall summary:\n")
                f.write(f"  Total input reads: {total_input}\n")
                if total_input > 0:
                    f.write(f"  Total passed reads: {total_passed} ({total_passed/total_input*100:.2f}%)\n")
                    f.write(f"  Total failed reads: {total_failed} ({total_failed/total_input*100:.2f}%)\n")
                else:
                    f.write(f"  Total passed reads: 0\n")
                    f.write(f"  Total failed reads: 0\n")
                f.write(f"\n")
                f.write(f"FINAL OUTPUT FILES:\n")
                f.write(f"-" * 60 + "\n")
                f.write(f"Standard files (platform-ordered: HiFi→ONT):\n")
                f.write(f"  {self.left_output}\n")
                f.write(f"  {self.right_output}\n")
                f.write(f"\n")
                f.write(f"Global merged files (platform-grouped):...\n")
                f.write(f"  {self.global_left_output}\n")
                f.write(f"    Structure: HiFi left → RC(HiFi right) → ONT left → RC(ONT right)\n")
                f.write(f"  {self.global_right_output}\n")
                f.write(f"    Structure: HiFi right → RC(HiFi left) → ONT right → RC(ONT left)\n")
                f.write(f"  Rationale: Platform grouping prioritizes high-quality HiFi reads\n")
                f.write(f"\n")
                f.write(f"Platform-specific files:\n")
                f.write(f"  {self.left_hifi_output}\n")
                f.write(f"  {self.right_hifi_output}\n")
                f.write(f"  {self.left_ont_output}\n")
                f.write(f"  {self.right_ont_output}\n")
                f.write(f"\n")
                f.write(f"Note: All temporary files (temp1.*, temp2.*) have been cleaned up.\n")
        
        print(f"[Part1-NEW] Log written to: {self.log_file}")
    
    def _print_summary(self):
        """Print summary statistics"""
        total = self.stats['total_reads']
        left = self.stats['left_telo']
        right = self.stats['right_telo']
        both = self.stats['both_telo']
        non = self.stats['non_telo']
        
        print(f"\n{'='*60}")
        print(f"SUMMARY (Window-based Algorithm)")
        print(f"{'='*60}")
        print(f"\nFirst-level filtering (TRC-based, window scanning):")
        print(f"  Total reads processed:           {total:>10}")
        
        if total > 0:
            print(f"  Left telomeric reads:            {left:>10} ({left/total*100:>5.1f}%)")
            print(f"  Right telomeric reads:           {right:>10} ({right/total*100:>5.1f}%)")
            print(f"  Both ends telomeric (excluded):  {both:>10} ({both/total*100:>5.1f}%)")
            print(f"  Non-telomeric reads:             {non:>10} ({non/total*100:>5.1f}%)")
        else:
            print(f"  Left telomeric reads:            {left:>10}")
            print(f"  Right telomeric reads:           {right:>10}")
            print(f"  Both ends telomeric (excluded):  {both:>10}")
            print(f"  Non-telomeric reads:             {non:>10}")
            print(f"\n  WARNING: No reads were processed. Check for errors above.")
        
        # Second filter summary
        if self.enable_second_filter:
            left_total = self.filter_stats['left']['total']
            left_passed = self.filter_stats['left']['passed']
            right_total = self.filter_stats['right']['total']
            right_passed = self.filter_stats['right']['passed']
            total_filter_input = left_total + right_total
            total_filter_passed = left_passed + right_passed
            
            print(f"\nSecond-level filtering (Motif count >= {self.min_motif_count}):")
            print(f"  Total input reads:               {total_filter_input:>10}")
            if total_filter_input > 0:
                print(f"  Total passed reads:              {total_filter_passed:>10} ({total_filter_passed/total_filter_input*100:>5.1f}%)")
                print(f"    Left passed:                   {left_passed:>10} ({left_passed/left_total*100:>5.1f}%)")
                print(f"    Right passed:                  {right_passed:>10} ({right_passed/right_total*100:>5.1f}%)")
            else:
                print(f"  Total passed reads:              {total_filter_passed:>10}")
        
        print(f"{'='*60}")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="TelSeekerPart1_new: Extract telomeric reads using window-based scanning",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        '--motif',
        required=True,
        help='Telomere repeat motif (e.g., TTAGGG)'
    )
    
    parser.add_argument(
        '--out',
        required=True,
        help='Output directory (where DEGAP placed reads_part/)'
    )
    
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=20,
        help='Number of threads for parallel processing (default: 20)'
    )
    
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.7,
        help='TRC threshold for telomeric classification (default: 0.7)'
    )
    
    parser.add_argument(
        '--check-length',
        type=int,
        default=1000,
        help='Length to check at each end in bp (default: 1000). Window flank = check_length / 2'
    )
    
    parser.add_argument(
        '--batch-size',
        type=int,
        default=1000,
        help='Number of reads to read at once for each batch (default: 1000)'
    )
    
    parser.add_argument(
        '--enable-second-filter',
        action='store_true',
        default=True,
        help='Enable second-level motif count filtering (default: enabled)'
    )
    
    parser.add_argument(
        '--disable-second-filter',
        action='store_true',
        help='Disable second-level motif count filtering'
    )
    
    parser.add_argument(
        '--min-motif-count',
        type=int,
        default=10,
        help='Minimum motif count for second-level filter (default: 10)'
    )
    
    parser.add_argument(
        '--overlapping',
        action='store_true',
        help='Use overlapping count in second-level filter (default: non-overlapping)'
    )
    
    args = parser.parse_args()
    
    # Handle filter enable/disable logic
    enable_second_filter = args.enable_second_filter and not args.disable_second_filter
    
    # Validate output directory
    if not os.path.exists(args.out):
        print(f"Error: Output directory not found: {args.out}")
        sys.exit(1)
    
    # Create extractor and run
    extractor = TeloReadsExtractor(
        output_dir=args.out,
        motif=args.motif,
        threads=args.threads,
        trc_threshold=args.threshold,
        check_length=args.check_length,
        batch_size=args.batch_size,
        enable_second_filter=enable_second_filter,
        min_motif_count=args.min_motif_count,
        overlapping=args.overlapping
    )
    
    extractor.run()


if __name__ == '__main__':
    main()
