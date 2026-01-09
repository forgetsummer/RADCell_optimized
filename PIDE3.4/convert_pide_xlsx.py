#!/usr/bin/env python3
"""
Convert PIDE3.4.xlsx to CSV format for C++ processing.

This script converts the Excel database file to a simple CSV format
that can be easily read by the C++ PIDEDataReader class.

Usage:
    python3 convert_pide_xlsx.py

Output:
    PIDE3.4.csv - CSV version of the main database
"""

import os
import sys

def main():
    try:
        import pandas as pd
    except ImportError:
        print("Error: pandas is required. Install with: pip3 install pandas openpyxl")
        sys.exit(1)
    
    try:
        import openpyxl
    except ImportError:
        print("Error: openpyxl is required. Install with: pip3 install openpyxl")
        sys.exit(1)
    
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    xlsx_file = os.path.join(script_dir, "PIDE3.4.xlsx")
    csv_file = os.path.join(script_dir, "PIDE3.4.csv")
    
    if not os.path.exists(xlsx_file):
        print(f"Error: {xlsx_file} not found!")
        sys.exit(1)
    
    print(f"Reading {xlsx_file}...")
    df = pd.read_excel(xlsx_file, engine='openpyxl')
    
    # Rename columns to match expected format (remove special characters)
    column_mapping = {
        '#ExpID': 'ExpID',
        '#Publication': 'PublicationID',
        '#IonExp': 'IonExpNum',
        '#PhotonExp': 'PhotonExpNum'
    }
    df = df.rename(columns=column_mapping)
    
    print(f"Loaded {len(df)} experiments")
    print(f"Columns: {list(df.columns)}")
    
    # Save to CSV
    print(f"Writing to {csv_file}...")
    df.to_csv(csv_file, index=False)
    
    print(f"Successfully created {csv_file}")
    print(f"\nDatabase summary:")
    print(f"  Total experiments: {len(df)}")
    print(f"  Unique cell lines: {df['Cells'].nunique()}")
    print(f"  Unique ions: {df['Ion'].nunique()}")
    
    # Show some statistics
    if 'CellClass' in df.columns:
        tumor = (df['CellClass'] == 't').sum()
        normal = (df['CellClass'] == 'n').sum()
        print(f"  Tumor cell experiments: {tumor}")
        print(f"  Normal cell experiments: {normal}")
    
    if 'CellOrigin' in df.columns:
        human = (df['CellOrigin'] == 'h').sum()
        rodent = (df['CellOrigin'] == 'r').sum()
        print(f"  Human cell experiments: {human}")
        print(f"  Rodent cell experiments: {rodent}")
    
    if 'LET' in df.columns:
        let_valid = df['LET'].dropna()
        print(f"  LET range: {let_valid.min():.2f} - {let_valid.max():.2f} keV/μm")
    
    # Print some example cell lines
    print(f"\nSample cell lines: {list(df['Cells'].unique()[:10])}")
    print(f"Sample ions: {list(df['Ion'].unique()[:10])}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
