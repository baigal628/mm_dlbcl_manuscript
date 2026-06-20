#!/usr/bin/env python3
"""
Calculate somatic hypermutation (SHM) rates for CDR1/CDR2 vs other V gene regions
from AIRR-formatted TSV files with IMGT-gapped alignments.
"""

import pandas as pd
import numpy as np
from scipy import stats
from typing import Dict, Tuple
import argparse
import sys
import re

dna_to_aa = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

def parse_imgt_positions(row: pd.Series) -> Dict[str, Tuple[int, int]]:
    """
    Parse IMGT numbering positions for different regions.
    IMGT positions (1-indexed, in gapped alignment):
    - FR1: 1-26
    - CDR1: 27-38
    - FR2: 39-55
    - CDR2: 56-65
    - FR3: 66-104
    - CDR3: 105+ (typically handled separately in AIRR format)
    The coordinate is based on amino acid. So need to convert to genomic coordinate.
    
    Returns positions as 0-indexed for Python string slicing.
    """
    regions = {
        'FR1': (0, 26*3),      # Positions 1-26 (0-indexed: 0-26, right open)
        'CDR1': (26*3, 38*3),    # Positions 27-38 (0-indexed: 26-38)
        'FR2': (38*3, 55*3),     # Positions 39-55 (0-indexed: 38-55)
        'CDR2': (55*3, 65*3),    # Positions 56-65 (0-indexed: 55-65)
        'FR3': (65*3, 104*3),    # Positions 66-104 (0-indexed: 65-104)
    }
    return regions

def count_mutations(seq_aligned: str, germline_aligned: str, 
                   start: int, end: int) -> Tuple[int, int]:
    """
    Count mutations and valid positions in a given region.
    
    Args:
        seq_aligned: IMGT-gapped sequence alignment
        germline_aligned: IMGT-gapped germline alignment
        start: Start position (0-indexed)
        end: End position (0-indexed, exclusive)
    
    Returns:
        (mutations, synonymous mutations, valid_positions): Number of mutations and synonmymous mutations, valid comparable positions
    """
    if pd.isna(seq_aligned) or pd.isna(germline_aligned):
        return np.nan, np.nan
    
    # Ensure strings are the same length
    min_len = min(len(seq_aligned), len(germline_aligned), end)
    if min_len <= start:
        return np.nan, np.nan
    
    end = min(end, min_len)
    
    mutations = 0
    synonymous_mutations = 0
    valid_positions = 0
    
    # concatenate the seq and germline. Only remove the positions with "."
    seq = seq_aligned[start:end].replace('.', '')
    germline = germline_aligned[start:end].replace('.', '')

    # Iterate through the region in codon triplets
    for i in range(0, len(seq), 3):
        # Skip codon contains indels
        codonWithIndel = False
        for j in range(i, i+3):
            if j >= len(seq): 
                break
            if seq[j] == '-' or germline[j] == '-':
                codonWithIndel = True
                break
        if codonWithIndel:
            continue

        for j in range(i, i+3):
            if j >= len(seq): 
                break
            seq_base = seq[j].upper()
            germ_base = germline[j].upper()
            
            # Skip ambiguous bases
            if seq_base in ['N', 'X'] or germ_base in ['N', 'X']:
                continue
            
            # Valid position for comparison
            valid_positions += 1
            
            # Count mutation if bases differ
            if seq_base != germ_base:
                # Check if it's a synonymous mutation
                seq_codon = seq[i:i+3].upper()
                germ_codon = germline[i:i+3].upper()
                if len(seq_codon) == 3 and len(germ_codon) == 3:
                    seq_aa = dna_to_aa.get(seq_codon, 'X')
                    germ_aa = dna_to_aa.get(germ_codon, 'X')
                    if seq_aa == germ_aa:
                        synonymous_mutations += 1
                mutations += 1
        
    return mutations, synonymous_mutations, valid_positions

def Lossos_test(sequence_alignment, germline_alignment, FR_mutations, FR_synonymous_mutations, CDR_mutations, CDR_synonymous_mutations, regions):
    # Calculate the baseline mutation probability that a random mutation replace the amino acid, Rf_FR, Rf_CDR
    Rf_FR = 0
    Rf_CDR = 0
    Rf_FR_total = 0 # total number of possible mutations in FR regions
    Rf_CDR_total = 0 # total number of possible mutations in CDR regions
    
    L_FR = regions['FR1'][1] - regions['FR1'][0] + regions['FR2'][1] - regions['FR2'][0] + regions['FR3'][1] - regions['FR3'][0]
    L_CDR = regions['CDR1'][1] - regions['CDR1'][0] + regions['CDR2'][1] - regions['CDR2'][0]
    
    tmp = L_FR+L_CDR
    L_FR = L_FR / tmp # relative fraction
    L_CDR = L_CDR / tmp

    for region_name, (start, end) in regions.items():
        germline = germline_alignment[start:end].replace('.', '')
        for i in range(0, start-end+1, 3):
            germ_codon = germ[i:i+3].upper()
            if len(germ_codon) < 3:
                continue
            if ("-" in germ_codon):
                continue
            for j in range(i, i+3):
                if j >= len(germ):
                    break
                germ_base = germ[j].upper()
                
                for base in ['A', 'C', 'G', 'T']:
                    if base == germ_base:
                        continue
                    mutated_codon = list(germ_codon)
                    mutated_codon[j-i] = base
                    mutated_codon = ''.join(mutated_codon)
                    germ_aa = dna_to_aa.get(germ_codon, 'X')
                    mutated_aa = dna_to_aa.get(mutated_codon, 'X')
                    if germ_aa != mutated_aa:
                         if region_name.startswith('FR'):
                            Rf_FR += 1
                         elif region_name.startswith('CDR'):
                            Rf_CDR += 1
                    if region_name.startswith('FR'):
                        Rf_FR_total += 1
                    elif region_name.startswith('CDR'):
                        Rf_CDR_total += 1 
    Rf_FR = Rf_FR / Rf_FR_total if Rf_FR_total > 0 else 0
    Rf_CDR = Rf_CDR / Rf_CDR_total if Rf_CDR_total > 0 else 0

    # core Lossos test
    # return two numbers: the p-value of CDR to have many replacement changes, and the p-value of FR to have very few replacement changes.
    def Test(s1, s2, r1, r2, p1, p2, q1, q2):
        ret = [0, 0]
        p = 0
        n = s1 + s2 + r1 + r2
        # Tet for CDR to have many synonymous mutation
        for k in range(r2, n+1):
            # enumerate all the combation for the variation S1,S2,R1, k where k+S1+S2+R2=n
            for S1 in range(n - k + 1):
                if (n - k - S1 <= 0):
                    continue
                for S2 in range(n - k - S1 + 1):                                                         
                    R1 = n - k - S1 - S2
                    if (R1 < 0):
                        continue
                    # calculate the multinomal of (n, R1, S1, R2, S2) with the probability of (p1, q1, p2, q2) 
                    x = stats.multinomial.pmf([R1, S1, k, S2], n, [p1, q1, p1, q2])
                    if (k == r2):
                        p += x/2
                    else:
                        p += x
        ret[0] = p
        
        p = 0
        # Test for FR to have very few synonymou mutation
        for k in range(0, r1+1):
            # enumerate all the combation for the variation k, S1,S2,R2 where k+S1+S2+R2=n
            for S1 in range(n - k + 1):
                if (n - k - S1 <= 0):
                    continue
                for S2 in range(n - k - S1 + 1):
                    R2 = n - k - S1 - S2
                    if (R2 < 0):
                        continue
                    # calculate the multinomal of (n, k, S1, S2, R2) with the probability of (p1, q1, p2, q2) 
                    x = stats.multinomial.pmf([k, S1, R2, S2], n, [p1, q1, p1, q2])
                    if (k == r1):
                        p += x/2
                    else:
                        p += x
        ret[1] = p
        return ret

    return Test(FR_synonymous_mutations, CDR_synonymous_mutations, FR_mutations - FR_synonymous_mutations, CDR_mutations - CDR_synonymous_mutations, Rf_FR * L_FR, Rf_CDR * L_CDR, (1-Rf_FR) * L_FR, (1-Rf_CDR) * L_CDR)

def ParseCigar(cigar):
    cigarFields = re.findall("\d+\w", cigar)
    ret = []
    for f in cigarFields:
        ret.append( (int(f[0:-1]), f[-1]) )
    return ret

def calculate_shm_rates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate SHM rates for each sequence.
    
    Args:
        df: DataFrame with AIRR format data
    
    Returns:
        DataFrame with added SHM rate columns
    """
    results = []
    
    for idx, row in df.iterrows():
        # Only check the IGH with complete VDJ
        if (row['complete_vdj'] == 'F'):
            continue 
        vfields = ParseCigar(row["v_cigar"])
        if ((len(vfields) >= 2 and vfields[0][1] == 'S' and vfields[1][1] == 'N')
           or vfields[0][1] == 'N'):
            continue
        jfields = ParseCigar(row["j_cigar"])
        if ((len(jfields) >= 2 and jfields[-1][1] == 'S' and jfields[-2][1] == 'N')
           or jfields[-1][1] == 'N'):
            continue

        if (not row['v_call'].startswith('IGH')
                or not row['j_call'].startswith('IGH')
                or not row['c_call'].startswith('IGH')):
            continue
        regions = parse_imgt_positions(row)
        
        # Get aligned sequences
        seq_aligned = row.get('sequence_alignment', '')
        germline_aligned = row.get('germline_alignment', '')
        
        # Calculate mutations for each region
        region_stats = {}
        for region_name, (start, end) in regions.items():
            mutations, synonymous_mutations, valid_pos = count_mutations(
                seq_aligned, germline_aligned, start, end
            )
            region_stats[region_name] = {
                'mutations': mutations,
                'synonymous_mutations': synonymous_mutations,
                'positions': valid_pos,
                'rate': mutations / valid_pos if valid_pos > 0 else np.nan
            }
        
        # Combine CDR regions
        cdr_mutations = (region_stats['CDR1']['mutations'] + 
                        region_stats['CDR2']['mutations'])
        cdr_synonymous = (region_stats['CDR1']['synonymous_mutations'] +
                        region_stats['CDR2']['synonymous_mutations'])
        cdr_positions = (region_stats['CDR1']['positions'] + 
                        region_stats['CDR2']['positions'])
        cdr_rate = cdr_mutations / cdr_positions if cdr_positions > 0 else np.nan
        cdr_syn_rate = cdr_synonymous / cdr_mutations if cdr_mutations > 0 else np.nan
        
        # Combine FR regions (FR1, FR2, FR3)
        fr_mutations = (region_stats['FR1']['mutations'] + 
                       region_stats['FR2']['mutations'] + 
                       region_stats['FR3']['mutations'])
        fr_synonymous = (region_stats['FR1']['synonymous_mutations'] + 
                       region_stats['FR2']['synonymous_mutations'] + 
                       region_stats['FR3']['synonymous_mutations'])
        fr_positions = (region_stats['FR1']['positions'] + 
                       region_stats['FR2']['positions'] + 
                       region_stats['FR3']['positions'])
        fr_rate = fr_mutations / fr_positions if fr_positions > 0 else np.nan
        fr_syn_rate = fr_synonymous / fr_mutations if fr_mutations > 0 else np.nan
        
        # Overall V gene (excluding CDR3)
        total_mutations = sum(s['mutations'] for s in region_stats.values())
        total_synonymous = sum(s['synonymous_mutations'] for s in region_stats.values())
        total_positions = sum(s['positions'] for s in region_stats.values())
        total_rate = total_mutations / total_positions if total_positions > 0 else np.nan
        total_syn_rate = total_synonymous / total_mutations if total_mutations > 0 else np.nan
       
        lossos_pvalues = Lossos_test(seq_aligned, germline_aligned, fr_mutations, fr_synonymous, cdr_mutations, cdr_synonymous, regions)

        # Store results
        result = {
            'sequence_id': row.get('sequence_id', idx),
            
            # Individual regions
            'FR1_mutations': region_stats['FR1']['mutations'],
            'FR1_synonymous': region_stats['FR1']['synonymous_mutations'],
            'FR1_positions': region_stats['FR1']['positions'],
            'FR1_rate': region_stats['FR1']['rate'],
            
            'CDR1_mutations': region_stats['CDR1']['mutations'],
            'CDR1_synonymous': region_stats['CDR1']['synonymous_mutations'],
            'CDR1_positions': region_stats['CDR1']['positions'],
            'CDR1_rate': region_stats['CDR1']['rate'],
            
            'FR2_mutations': region_stats['FR2']['mutations'],
            'FR2_synonymous': region_stats['FR2']['synonymous_mutations'],
            'FR2_positions': region_stats['FR2']['positions'],
            'FR2_rate': region_stats['FR2']['rate'],
            
            'CDR2_mutations': region_stats['CDR2']['mutations'],
            'CDR2_synonymous': region_stats['CDR2']['synonymous_mutations'],
            'CDR2_positions': region_stats['CDR2']['positions'],
            'CDR2_rate': region_stats['CDR2']['rate'],
            
            'FR3_mutations': region_stats['FR3']['mutations'],
            'FR3_synonymous': region_stats['FR3']['synonymous_mutations'],
            'FR3_positions': region_stats['FR3']['positions'],
            'FR3_rate': region_stats['FR3']['rate'],
            
            # Combined regions
            'CDR_mutations': cdr_mutations,
            'CDR_synonymous': cdr_synonymous,
            'CDR_positions': cdr_positions,
            'CDR_rate': cdr_rate,
            'CDR_syn_rate': cdr_syn_rate,
            
            'FR_mutations': fr_mutations,
            'FR_synonymous': fr_synonymous,
            'FR_positions': fr_positions,
            'FR_rate': fr_rate,
            'FR_syn_rate': fr_syn_rate,
            
            'total_V_mutations': total_mutations,
            'total_V_synonymous': total_synonymous,
            'total_V_positions': total_positions,
            'total_V_rate': total_rate,
            'total_V_syn_rate': total_syn_rate,

            'CDR_Lossos_pvalue': lossos_pvalues[0],
            'FR_Lossos_pvalue': lossos_pvalues[1],
        }
        
        results.append(result)
    print(f"Valid {len(results)} items") 
    return pd.DataFrame(results)

def calculate_summary_statistics(shm_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate summary statistics across all sequences.
    """
    summary = {
        'Region': ['CDR1', 'CDR2', 'CDR1+CDR2', 'CDR1+CDR2 synonymous', 'FR1', 'FR2', 'FR3', 'FR1+FR2+FR3', 'FR1+FR2+FR3 synonymous', 'Total V', 'Total V synonymous', 'Lossos CDR p-value', 'Lossos FR p-value'],
        'Mean_Rate': [
            shm_df['CDR1_rate'].mean(),
            shm_df['CDR2_rate'].mean(),
            shm_df['CDR_rate'].mean(),
            shm_df['CDR_syn_rate'].mean(),
            shm_df['FR1_rate'].mean(),
            shm_df['FR2_rate'].mean(),
            shm_df['FR3_rate'].mean(),
            shm_df['FR_rate'].mean(),
            shm_df['FR_syn_rate'].mean(),
            shm_df['total_V_rate'].mean(),
            shm_df['total_V_syn_rate'].mean(),
            shm_df['CDR_Lossos_pvalue'].mean(),
            shm_df['FR_Lossos_pvalue'].mean(),
        ],
        'Median_Rate': [
            shm_df['CDR1_rate'].median(),
            shm_df['CDR2_rate'].median(),
            shm_df['CDR_rate'].median(),
            shm_df['CDR_syn_rate'].median(),
            shm_df['FR1_rate'].median(),
            shm_df['FR2_rate'].median(),
            shm_df['FR3_rate'].median(),
            shm_df['FR_rate'].median(),
            shm_df['FR_syn_rate'].median(),
            shm_df['total_V_rate'].median(),
            shm_df['total_V_syn_rate'].median(),
            shm_df['CDR_Lossos_pvalue'].median(),
            shm_df['FR_Lossos_pvalue'].median(),
        ],
        'Std_Rate': [
            shm_df['CDR1_rate'].std(),
            shm_df['CDR2_rate'].std(),
            shm_df['CDR_rate'].std(),
            shm_df['CDR_syn_rate'].std(),
            shm_df['FR1_rate'].std(),
            shm_df['FR2_rate'].std(),
            shm_df['FR3_rate'].std(),
            shm_df['FR_rate'].std(),
            shm_df['FR_syn_rate'].std(),
            shm_df['total_V_rate'].std(),
            shm_df['total_V_syn_rate'].std(),
            shm_df['CDR_Lossos_pvalue'].std(),
            shm_df['FR_Lossos_pvalue'].std(),
        ],
        'Total_Mutations': [
            shm_df['CDR1_mutations'].sum(),
            shm_df['CDR2_mutations'].sum(),
            shm_df['CDR_mutations'].sum(),
            shm_df['CDR_synonymous'].sum(),
            shm_df['FR1_mutations'].sum(),
            shm_df['FR2_mutations'].sum(),
            shm_df['FR3_mutations'].sum(),
            shm_df['FR_mutations'].sum(),
            shm_df['FR_synonymous'].sum(),
            shm_df['total_V_mutations'].sum(),
            shm_df['total_V_synonymous'].sum(),
            np.nan, # Lossos p-values don't have "mutations"
            np.nan,
        ],
        'Total_Positions': [
            shm_df['CDR1_positions'].sum(),
            shm_df['CDR2_positions'].sum(),
            shm_df['CDR_positions'].sum(),
            shm_df['CDR_mutations'].sum(),
            shm_df['FR1_positions'].sum(),
            shm_df['FR2_positions'].sum(),
            shm_df['FR3_positions'].sum(),
            shm_df['FR_positions'].sum(),
            shm_df['FR_mutations'].sum(),  
            shm_df['total_V_positions'].sum(),
            shm_df['total_V_mutations'].sum(),  
            np.nan, # Lossos p-values don't have "positions"
            np.nan,
        ],
    }
    
    summary_df = pd.DataFrame(summary)
    summary_df['Overall_Rate'] = (summary_df['Total_Mutations'] / 
                                   summary_df['Total_Positions'])
    
    return summary_df

def main():
    parser = argparse.ArgumentParser(
        description='Calculate somatic hypermutation rates from AIRR-formatted TSV files'
    )
    parser.add_argument('input', help='Input TSV file in AIRR format')
    parser.add_argument('-o', '--output', default='shm_results.tsv',
                       help='Output file for per-sequence SHM rates (default: shm_results.tsv)')
    parser.add_argument('-s', '--summary', default='shm_summary.tsv',
                       help='Output file for summary statistics (default: shm_summary.tsv)')
    parser.add_argument('--seq-col', default='sequence_alignment',
                       help='Column name for sequence alignment (default: sequence_alignment)')
    parser.add_argument('--germ-col', default='germline_alignment',
                       help='Column name for germline alignment (default: germline_alignment)')
    
    args = parser.parse_args()
    
    # Read input file
    print(f"Reading input file: {args.input}")
    try:
        df = pd.read_csv(args.input, sep='\t')
        print(f"Loaded {len(df)} sequences")
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check required columns
    if args.seq_col not in df.columns or args.germ_col not in df.columns:
        print(f"Error: Required columns '{args.seq_col}' and/or '{args.germ_col}' not found", 
              file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)
    
    # Calculate SHM rates
    print("Calculating SHM rates...")
    shm_results = calculate_shm_rates(df)
    
    # Merge with original data
    output_df = shm_results # pd.concat([df, shm_results.drop('sequence_id', axis=1)], axis=1)
    
    # Save per-sequence results
    print(f"Writing per-sequence results to: {args.output}")
    output_df.to_csv(args.output, sep='\t', index=False)
    
    # Calculate and save summary statistics
    print("Calculating summary statistics...")
    summary_df = calculate_summary_statistics(shm_results)
    print(f"Writing summary statistics to: {args.summary}")
    summary_df.to_csv(args.summary, sep='\t', index=False)
    
    # Print summary to console
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("\n")
    print(f"CDR1+CDR2 vs FR1+FR2+FR3 ratio: "
          f"{summary_df.loc[2, 'Overall_Rate'] / summary_df.loc[6, 'Overall_Rate']:.3f}")
    print("="*80)
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
