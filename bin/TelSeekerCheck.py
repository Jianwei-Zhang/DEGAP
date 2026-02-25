#!/usr/bin/env python3
"""
TelSeekerCheck.py - Check if chromosome ends have reached telomeres

This script analyzes a genome assembly and identifies which chromosome ends
have successfully reached telomeric sequences and which need telomere extension.

Modified to use window-based scanning method:
- Extract 5kb from each chromosome end (default, to avoid ITS region interference)
- Use window markers (motif*2) to scan for telomeric regions
- Create windows (500bp upstream + marker + 500bp downstream)
- Calculate TRC for each window and take the maximum
- If max TRC >= 0.7, consider the end as telomeric

Usage:
    python TelSeekerCheck.py --genome genome.fa --motif TTAGGG --out output_dir

Output:
    output_dir/genome.telomere.check/
        ├── genome.telomere.check.csv          # Summary report
        ├── genome.telomere.check.left.2kb.fa  # Left end sequences (2kb)
        ├── genome.telomere.check.right.2kb.fa # Right end sequences (2kb)
        └── need_extension_chr_end.txt         # List of ends needing extension
"""

import os
import sys
import re
import argparse
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ============================================================================
# TRC Calculation Functions (Topsicle-based algorithm)
# ============================================================================

def _pattern_scramble_telo(pattern: str, cut_length: int) -> List[str]:
    """
    Generate all k-mer permutations of telomere motif (including reverse complement).
    
    Args:
        pattern: Telomere motif (e.g., 'TTAGGG')
        cut_length: Length of k-mers to generate
        
    Returns:
        Sorted list of unique k-mers
    """
    pattern = (pattern or "").upper()
    if not pattern:
        return []
    
    cut_length = max(1, cut_length)
    extended = pattern * 2
    kmers = set()
    
    # Generate cyclic rotations
    for i in range(len(pattern)):
        fragment = extended[i : i + cut_length]
        if len(fragment) == cut_length:
            kmers.add(fragment)
    
    # Add reverse complements
    kmers.update(str(Seq(kmer).reverse_complement()) for kmer in list(kmers))
    return sorted(kmers)


def patterns_to_search_topsicle(motif: str, cut_length: int) -> List[str]:
    """
    Generate regex patterns for Topsicle-based telomere detection.

    This non-overlapping implementation mirrors the original Topsicle behaviour:
    each pattern is matched in a classic left-to-right scan without lookahead,
    so successive hits do not share bases.

    Args:
        motif: Telomere motif
        cut_length: Length of k-mers

    Returns:
        List of escaped regex patterns for non-overlapping matching
    """
    kmers = _pattern_scramble_telo(motif, cut_length)
    return [re.escape(kmer) for kmer in kmers]


def calculate_trc(sequence: str, motif: str, kmer_length: int, check_length: int = 1000) -> float:
    """
    Calculate Telomere Repeat Content (TRC) using Topsicle method.

    This function implements the exact algorithm from Topsicle:
    1. Generate all telomere k-mer patterns (including reverse complements)
    2. Calculate match density for each pattern
    3. Return the maximum density as TRC score

    Args:
        sequence: DNA sequence to check
        motif: Telomeric motif (e.g., 'TTAGGG')
        kmer_length: Length of k-mers to use (typically motif_length - 2)
        check_length: Length of sequence to check from the start

    Returns:
        TRC score between 0 and infinity (typically 0-1 range)
    """
    # Truncate sequence to check_length
    seq_to_check = sequence[:min(len(sequence), check_length)].upper()

    if len(seq_to_check) < kmer_length:
        return 0.0

    # Generate Topsicle patterns (non-overlapping)
    pattern_all = patterns_to_search_topsicle(motif, kmer_length)
    compiled_patterns = [re.compile(pattern) for pattern in pattern_all]

    # Calculate TRC using Topsicle method: max density across all patterns
    # Ensure ratio_perfect_hit is at least 1.0
    ratio_perfect_hit = max(1.0, len(seq_to_check) / len(motif))
    max_trc = 0.0

    for pattern in compiled_patterns:
        matches = len(pattern.findall(seq_to_check))
        pattern_trc = matches / ratio_perfect_hit if ratio_perfect_hit > 0 else 0.0
        max_trc = max(max_trc, pattern_trc)

    return max_trc


def find_all_occurrences(seq: str, pattern: str) -> List[int]:
    """
    Find all occurrences of pattern in sequence (overlapping allowed for finding positions).

    Args:
        seq: DNA sequence to search
        pattern: Pattern to find

    Returns:
        List of starting positions
    """
    positions = []
    start = 0
    while True:
        pos = seq.find(pattern, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1  # Allow overlapping search
    return positions


def calculate_trc_window_based(sequence: str, motif: str, kmer_length: int,
                                window_flank: int = 500) -> Tuple[float, int]:
    """
    Calculate maximum TRC using window-based scanning method.

    This method:
    1. Generates window markers (motif*2 and its reverse complement)
    2. Finds all marker positions in the sequence
    3. Creates windows around each marker (flank + marker + flank)
    4. Calculates TRC for each window
    5. Returns the maximum TRC found

    Args:
        sequence: DNA sequence to check (e.g., 5kb from chromosome end)
        motif: Telomeric motif (e.g., 'TTAGGG')
        kmer_length: Length of k-mers for TRC calculation
        window_flank: Flank size on each side of marker (default: 500bp)

    Returns:
        Tuple of (max_trc, num_windows_found)
    """
    seq_str = sequence.upper()

    # Generate window markers (motif*2)
    pat_right = motif * 2  # e.g., TTAGGGTTAGGG
    pat_left = str(Seq(pat_right).reverse_complement())  # e.g., CCCTAACCCTAA

    marker_len = len(pat_right)
    target_window_size = 2 * window_flank + marker_len

    # Find all marker positions
    right_positions = find_all_occurrences(seq_str, pat_right)
    left_positions = find_all_occurrences(seq_str, pat_left)

    total_windows = len(right_positions) + len(left_positions)

    if total_windows == 0:
        return 0.0, 0

    max_trc = 0.0

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
            # Don't pass check_length - use the full window sequence
            window_trc = calculate_trc(window_seq, motif, kmer_length, len(window_seq))
            max_trc = max(max_trc, window_trc)

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
            # Don't pass check_length - use the full window sequence
            window_trc = calculate_trc(window_seq, motif, kmer_length, len(window_seq))
            max_trc = max(max_trc, window_trc)

    return max_trc, total_windows


# ============================================================================
# Telomere Checker Main Class
# ============================================================================


class TelomereChecker:
    """Check chromosome ends for telomeric content using window-based scanning"""

    def __init__(self, genome_file: str, motif: str, output_dir: str,
                 trc_threshold: float = 0.7, check_length: int = 5000,
                 extract_length: int = 2000, window_flank: int = 500):
        """
        Initialize TelomereChecker

        Args:
            genome_file: Path to input genome FASTA file
            motif: Telomere repeat motif (e.g., 'TTAGGG')
            output_dir: Output directory
            trc_threshold: TRC threshold to consider telomeric (default: 0.7)
            check_length: Length to extract from chromosome end for window scanning (default: 5000bp)
            extract_length: Length to extract for output sequences (default: 2000bp)
            window_flank: Flank size for window creation (default: 500bp)
        """
        self.genome_file = genome_file
        self.motif = motif.upper()
        self.output_dir = output_dir
        self.trc_threshold = trc_threshold
        self.check_length = check_length
        self.extract_length = extract_length
        self.window_flank = window_flank
        self.kmer_length = len(motif) - 2

        # Window markers (motif*2)
        self.pat_right = self.motif * 2
        self.pat_left = str(Seq(self.pat_right).reverse_complement())

        # Results storage
        self.results = []
        self.left_sequences = []
        self.right_sequences = []
        
    def run(self):
        """Run the complete telomere checking workflow"""
        print(f"[TelomereChecker] Starting telomere check (window-based method)...")
        print(f"  Genome: {self.genome_file}")
        print(f"  Motif: {self.motif}")
        print(f"  Window markers:")
        print(f"    Right (forward): {self.pat_right}")
        print(f"    Left (reverse):  {self.pat_left}")
        print(f"  TRC threshold: {self.trc_threshold}")
        print(f"  Check length: {self.check_length} bp (extract from each end)")
        print(f"  Window flank: {self.window_flank} bp")
        print(f"  Extract length: {self.extract_length} bp (for output)")
        print()

        # Step 1: Create output directory
        self._create_output_directory()

        # Step 2: Process genome
        self._process_genome()

        # Step 3: Write results
        self._write_csv_report()
        self._write_extension_list()
        self._write_left_sequences()
        self._write_right_sequences()

        # Step 4: Print summary
        self._print_summary()

        print(f"\n[TelomereChecker] Complete! Results saved to: {self.check_dir}")
        
    def _create_output_directory(self):
        """Create output directory structure"""
        # Fixed directory name: genome.telomere.check
        # This is hardcoded regardless of input genome file name
        self.check_dir = os.path.join(self.output_dir, "genome.telomere.check")
        os.makedirs(self.check_dir, exist_ok=True)
        
        print(f"[TelomereChecker] Output directory: {self.check_dir}")
        
    def _process_genome(self):
        """Process genome and check each chromosome end"""
        print(f"[TelomereChecker] Processing genome...")
        
        try:
            with open(self.genome_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    chr_name = record.id
                    seq_str = str(record.seq)
                    chr_length = len(seq_str)
                    
                    print(f"  Processing {chr_name} ({chr_length} bp)...")
                    
                    # Check left end
                    left_result = self._check_end(chr_name, seq_str, 'Left')
                    self.results.append(left_result)
                    
                    # Extract left sequence
                    left_seq = self._extract_left_sequence(chr_name, seq_str)
                    if left_seq:
                        self.left_sequences.append(left_seq)
                    
                    # Check right end
                    right_result = self._check_end(chr_name, seq_str, 'Right')
                    self.results.append(right_result)
                    
                    # Extract right sequence
                    right_seq = self._extract_right_sequence(chr_name, seq_str)
                    if right_seq:
                        self.right_sequences.append(right_seq)
                        
        except FileNotFoundError:
            print(f"[TelomereChecker] Error: Genome file not found: {self.genome_file}")
            sys.exit(1)
        except Exception as e:
            print(f"[TelomereChecker] Error processing genome: {e}")
            sys.exit(1)
            
    def _check_end(self, chr_name: str, seq_str: str, end: str) -> dict:
        """
        Check one chromosome end for telomeric content using window-based scanning.

        Args:
            chr_name: Chromosome name
            seq_str: Full chromosome sequence
            end: 'Left' or 'Right'

        Returns:
            dict with check results including max_trc and num_windows
        """
        # Extract sequence to check (5kb by default)
        if end == 'Left':
            check_seq = seq_str[:min(self.check_length, len(seq_str))]
        else:  # Right
            check_seq = seq_str[-min(self.check_length, len(seq_str)):]

        # Calculate TRC using window-based method
        max_trc, num_windows = calculate_trc_window_based(
            check_seq,
            self.motif,
            self.kmer_length,
            self.window_flank
        )

        # Determine status based on maximum TRC across all windows
        status = "telomeric" if max_trc >= self.trc_threshold else "untelomeric"

        return {
            'chr': chr_name,
            'end': end,
            'status': status,
            'trc': max_trc,
            'num_windows': num_windows,
            'check_length': len(check_seq)
        }
    
    def _extract_left_sequence(self, chr_name: str, seq_str: str) -> SeqRecord:
        """Extract left end sequence"""
        extract_len = min(self.extract_length, len(seq_str))
        left_seq = seq_str[:extract_len]
        
        return SeqRecord(
            seq=Seq(left_seq),
            id=f"{chr_name}_L",
            description=f"Left end of {chr_name}, {extract_len} bp"
        )
    
    def _extract_right_sequence(self, chr_name: str, seq_str: str) -> SeqRecord:
        """Extract right end sequence"""
        extract_len = min(self.extract_length, len(seq_str))
        right_seq = seq_str[-extract_len:]
        
        return SeqRecord(
            seq=Seq(right_seq),
            id=f"{chr_name}_R",
            description=f"Right end of {chr_name}, {extract_len} bp"
        )
    
    def _write_csv_report(self):
        """Write CSV report with window-based results"""
        # Fixed filename: genome.telomere.check.csv
        csv_file = os.path.join(self.check_dir, "genome.telomere.check.csv")

        print(f"\n[TelomereChecker] Writing CSV report: {csv_file}")

        with open(csv_file, 'w') as f:
            # Write results with additional window information
            for result in self.results:
                # Format: Chr01,Left,telomeric,TRC=0.98,Windows=5,CheckLen=5000
                line = (f"{result['chr']},{result['end']},{result['status']},"
                       f"TRC={result['trc']:.2f},"
                       f"Windows={result['num_windows']},"
                       f"CheckLen={result['check_length']}\n")
                f.write(line)

        print(f"  Written {len(self.results)} records")
    
    def _write_extension_list(self):
        """Write list of chromosome ends that need extension"""
        extension_file = os.path.join(self.check_dir, "need_extension_chr_end.txt")
        
        print(f"[TelomereChecker] Writing extension list: {extension_file}")
        
        # Collect untelomeric ends
        untelomeric_ends = []
        for result in self.results:
            if result['status'] == 'untelomeric':
                # Format: Chr01.L or Chr01.R
                end_suffix = 'L' if result['end'] == 'Left' else 'R'
                untelomeric_ends.append(f"{result['chr']}.{end_suffix}")
        
        # Write to file
        with open(extension_file, 'w') as f:
            for end in untelomeric_ends:
                f.write(f"{end}\n")
        
        print(f"  Written {len(untelomeric_ends)} ends requiring extension")
    
    def _write_left_sequences(self):
        """Write left end sequences"""
        # Fixed filename: genome.telomere.check.left.2kb.fa
        left_file = os.path.join(self.check_dir, "genome.telomere.check.left.2kb.fa")
        
        print(f"[TelomereChecker] Writing left sequences: {left_file}")
        
        with open(left_file, 'w') as f:
            SeqIO.write(self.left_sequences, f, 'fasta')
        
        print(f"  Written {len(self.left_sequences)} sequences")
    
    def _write_right_sequences(self):
        """Write right end sequences"""
        # Fixed filename: genome.telomere.check.right.2kb.fa
        right_file = os.path.join(self.check_dir, "genome.telomere.check.right.2kb.fa")
        
        print(f"[TelomereChecker] Writing right sequences: {right_file}")
        
        with open(right_file, 'w') as f:
            SeqIO.write(self.right_sequences, f, 'fasta')
        
        print(f"  Written {len(self.right_sequences)} sequences")
    
    def _print_summary(self):
        """Print summary statistics with window information"""
        total_ends = len(self.results)
        telomeric_count = sum(1 for r in self.results if r['status'] == 'telomeric')
        untelomeric_count = sum(1 for r in self.results if r['status'] == 'untelomeric')

        print(f"\n{'='*60}")
        print(f"SUMMARY (Window-based Method)")
        print(f"{'='*60}")
        print(f"Total chromosome ends checked: {total_ends}")
        print(f"  Telomeric ends:   {telomeric_count} ({telomeric_count/total_ends*100:.1f}%)")
        print(f"  Untelomeric ends: {untelomeric_count} ({untelomeric_count/total_ends*100:.1f}%)")
        print(f"{'='*60}")

        if untelomeric_count > 0:
            print(f"\nChromosome ends requiring telomere extension:")
            for result in self.results:
                if result['status'] == 'untelomeric':
                    print(f"  - {result['chr']} ({result['end']}): "
                          f"TRC={result['trc']:.2f}, Windows={result['num_windows']}")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Check if chromosome ends have reached telomeres (window-based method)",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example usage:
  python TelSeekerCheck.py --genome genome.fa --motif TTAGGG --out results

Window-based method:
  - Extracts 5kb from each chromosome end (default, to avoid ITS interference)
  - Uses window markers (motif*2) to scan for telomeric regions
  - Creates windows (500bp upstream + marker + 500bp downstream)
  - Calculates TRC for each window and takes the maximum
  - If max TRC >= 0.7, considers the end as telomeric

Output structure:
  results/genome.telomere.check/
    ├── genome.telomere.check.csv          # Summary report with window info
    ├── need_extension_chr_end.txt         # List of ends needing extension
    ├── genome.telomere.check.left.2kb.fa  # Left end sequences
    └── genome.telomere.check.right.2kb.fa # Right end sequences
        """
    )

    parser.add_argument(
        '--genome',
        required=True,
        help='Input genome FASTA file'
    )

    parser.add_argument(
        '--motif',
        default='TTAGGG',
        help='Telomere repeat motif (default: TTAGGG)'
    )

    parser.add_argument(
        '--out',
        required=True,
        help='Output directory'
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
        default=5000,
        help='Length to extract from chromosome end for window scanning in bp (default: 5000)'
    )

    parser.add_argument(
        '--extract-length',
        type=int,
        default=2000,
        help='Length to extract for output sequences in bp (default: 2000)'
    )

    parser.add_argument(
        '--window-flank',
        type=int,
        default=500,
        help='Flank size for window creation in bp (default: 500)'
    )

    args = parser.parse_args()

    # Validate genome file
    if not os.path.exists(args.genome):
        print(f"Error: Genome file not found: {args.genome}")
        sys.exit(1)

    # Create checker and run
    checker = TelomereChecker(
        genome_file=args.genome,
        motif=args.motif,
        output_dir=args.out,
        trc_threshold=args.threshold,
        check_length=args.check_length,
        extract_length=args.extract_length,
        window_flank=args.window_flank
    )

    checker.run()


if __name__ == '__main__':
    main()
