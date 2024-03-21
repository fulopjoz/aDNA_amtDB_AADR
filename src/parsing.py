#!/usr/bin/env python3
# Author: Jozef Fulop
# Institution: UCT in Prague

import pandas as pd
import re

def parse_publication_name(reference):
    # Tries to extract the first capitalized name and year to format the publication name.
    match_first_name = re.search(r"[A-Z][a-z]+", reference)
    match_year = re.search(r"\d{4}", reference)
    if match_first_name and match_year:
        first_name = match_first_name.group()
        year = match_year.group()
        return f"{first_name} et al. {year}"
    return ""

def match_reich_metadata(reich_meta_file, missing_ids, mitopatho_csv):
    print("Reading Reich metadata...")
    reich_meta = pd.read_csv(reich_meta_file, sep='\t', header=0, low_memory=False)
    mitopatho = pd.read_csv(mitopatho_csv)
    print(f"Reich metadata has {len(reich_meta)} rows.")
    print(f"MitoPatho data has {len(mitopatho)} rows.")
    
    missing_ids_set = set(missing_ids)
    matched_data = []
    
    for _, row in reich_meta.iterrows():
        master_id = row['Master ID']
        if master_id not in missing_ids_set:
            continue
        
        # Attempt to parse the calibrated radiocarbon age and archaeological context range
        date_detail = row.get('Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE', '')
        date_range_match = re.search(r'(\d{4})\u00b1(\d{2})', date_detail)
        c14_lab_code_match = re.search(r'Ua-(\w+)', date_detail)
        year_from, year_to, c14_lab_code = None, None, 'Unknown'

        if date_range_match:
            bp, error = date_range_match.groups()
            year_from = -int(bp)  # Assume the presence of 'BP' indicates Before Present
            year_to = -(int(bp) - int(error) * 2)
        if c14_lab_code_match:
            c14_lab_code = c14_lab_code_match.group(1)

        c14_sample_tag = 1 if year_to and year_to > 1950 else 0

        # Extract data from MitoPatho if available
        mito_filtered = mitopatho[mitopatho['ID'] == master_id]
        mito_agg = {col: ";".join(mito_filtered[col].fillna('<NA>').astype(str).unique()) for col in ['Polymorphism', 'Position', 'Locus', 'Diseases', 'Status', 'Homoplasmy', 'Heteroplasmy']}
        
        matched_data.append({
            'identifier':               row['Master ID'], 
            'alternative_identifiers':  row.get('Genetic ID', ""), 
            'country':                  row.get('Political Entity', ""), 
            'continent':                "",
            'region':                   "",
            'culture':                  row.get('Group ID', ""),
            'epoch':                    "",
            'group':                    row.get('Group ID', ""),
            'comment':                  "",
            'latitude':                 row.get('Lat.', ""),
            'longitude':                row.get('Long.', ""), 
            'sex':                      row.get('Molecular Sex', ""), 
            'site':                     row.get('Locality', ""), 
            'site_detail':              "",
            'mt_hg':                    row.get('mtDNA haplogroup if >2x or published', ""),
            'ychr_hg':                  row.get('Y haplogroup (manual curation in ISOGG format)', ""),
            'year_from':                year_from, 
            'year_to':                  year_to, 
            'date_detail':              date_detail,
            'bp':                       bp if date_range_match else '', 
            'c14_lab_code':             c14_lab_code,
            'reference_name':           parse_publication_name(row.get('Publication', "")),
            'c14_sample_tag':           c14_sample_tag,
            'ychr_snps':                row.get('Y haplogroup (manual curation in terminal mutation format)', ""),
            'avg_coverage':             row.get('mtDNA coverage (merged data)', ""),
            'sequence_source':          'fasta',
            'mitopatho_alleles':        mito_agg.get('Polymorphism', ''),
            'mitopatho_positions':      mito_agg.get('Position', ''),
            'mitopatho_locus':          mito_agg.get('Locus', ''),
            'mitopatho_diseases':       mito_agg.get('Diseases', ''),
            'mitopatho_statuses':       mito_agg.get('Status', ''),
            'mitopatho_homoplasms':     mito_agg.get('Homoplasmy', ''),
            'mitopatho_heteroplasms':   mito_agg.get('Heteroplasmy', '')
        })

    print(f"Processed {len(matched_data)} out of {len(missing_ids)} missing IDs.")
    matched_df = pd.DataFrame(matched_data)
    return matched_df


