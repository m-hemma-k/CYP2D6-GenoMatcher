#!/usr/bin/env python3

import io
import re
import pandas as pd  # type: ignore

def read_vcf(path):
    """
    Reads a VCF file and returns its contents as a Pandas DataFrame.
    Ignores metadata lines starting with '##'.

    Args:
        path (str): Path to the VCF file.

    Returns:
        pd.DataFrame: A DataFrame containing the VCF file data.
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def extract_rs_genotypes(diplotypes_df, vcf_df):
    """
    Extract rsIDs from the diplotypes DataFrame and their corresponding genotypes 
    from the VCF DataFrame.

    Args:
        diplotypes_df (pd.DataFrame): DataFrame containing rsIDs as column headers.
        vcf_df (pd.DataFrame): VCF DataFrame containing variant information.

    Returns:
        dict: A dictionary with rsIDs as keys and genotypes as values.
    """
    genotypes = {}

    # Step 1: Extract rsIDs from diplotypes_df
    rs_pattern = re.compile(r'^rs\d+$')
    valid_rsids = [col for col in diplotypes_df.columns if rs_pattern.match(col)]

    # Step 2: Extract genotypes from VCF DataFrame only for valid rsIDs
    for _, row in vcf_df.iterrows():
        if pd.notna(row['ID']) and row['ID'] in valid_rsids:  # Check if rsID is valid
            try:
                # Locate GT (genotype) field position in FORMAT
                format_fields = row['FORMAT'].split(':')
                gt_index = format_fields.index('GT')  # Find position of GT
                
                # Extract genotype from the first sample column (10th column)
                sample_data = row.iloc[9]  # Assuming first sample column
                sample_fields = sample_data.split(':')
                gt_value = sample_fields[gt_index]
                
                # Convert GT (e.g., '0/1') to allele letters
                alleles = [row['REF']] + row['ALT'].split(',')
                genotype = '/'.join([alleles[int(i)] for i in gt_value.replace('|', '/').split('/') if i.isdigit()])
                
                # Add to dictionary
                genotypes[row['ID']] = genotype
            except (IndexError, ValueError) as e:
                print(f"Error processing row {row['ID']}: {e}")
                continue
    
    return genotypes


# Funktion für den Abgleich und Sortierung nach Ranking
def match_and_sort(diplotypes_df, sample_data):
    # Entferne die Spalte 'Ranking' für den Abgleich
    diplotypes_df_filtered = diplotypes_df.drop(columns=['Ranking'])

    # Abgleich der Zeilen mit sample_data
    matching_rows = diplotypes_df_filtered[
        diplotypes_df_filtered.apply(
            lambda row: all(row[key] == value for key, value in sample_data.items() if key in row.index), axis=1
        )
    ]

    # Füge die Ranking-Spalte hinzu und sortiere die Ergebnisse
    matching_rows['Ranking'] = diplotypes_df.loc[matching_rows.index, 'Ranking']
    matching_rows = matching_rows.sort_values(by='Ranking', ascending=False)

    # Nur die Genotype-Spalte ausgeben
    return matching_rows['Genotype']


def print_matches(diplotypes_df, matches):
    """
    Prints matches sorted by 'Ranking' value.

    Args:
        diplotypes_df (pd.DataFrame): Reference genotype DataFrame.
        matches (list): List of matching row indices.
    """
    matches_sorted = sorted(matches, key=lambda x: diplotypes_df.loc[x, 'Ranking'], reverse=True)
    print("Matches sorted by probability (Ranking) descending:")
    for match in matches_sorted:
        tier_value = diplotypes_df.loc[match, 'Ranking']
        cnv_value = diplotypes_df.loc[match, 'CNV']
        print(f"100% Match: Row {match}, Ranking: {tier_value}, CNV: {cnv_value}")


def main():
    """
    Main function to load VCF data, extract genotypes, compare against reference,
    and print matching rows.
    """
    # Load the reference DataFrame
    file_path = '../00-preprocessing_data/output/CYP2D6_2024-12-16.pkl'
    diplotypes_df = pd.read_pickle(file_path)

    # Load the VCF file
    sample_filepath = './input/Test_Input_1_1.vcf'
    vcf_data = read_vcf(sample_filepath)
    # print(vcf_data)
    # Save as CSV
    # diplotypes_csv_path = './N8A1499_PharmCatInput.csv'
    # diplotypes_df.to_csv(diplotypes_csv_path, index=False)

    sample_data = extract_rs_genotypes(diplotypes_df, vcf_data)

    # Extrahiere CNV-Wert aus vcf_data
    CNV = vcf_data.loc[vcf_data['ID'] == 'CYP2D6_CNV', 'PharmCAT']
    CNV_value = CNV.iloc[0] if not CNV.empty else None
    sample_data = extract_rs_genotypes(diplotypes_df, vcf_data)
    sample_data['CNV'] = CNV_value

    # Ausgabe des Ergebnisses
    print(sample_data)

    # Entferne die 'Ranking'-Spalte aus diplotypes_df
    diplotypes_df_filtered = diplotypes_df.drop(columns=['Ranking'], errors='ignore')

    # Abgleich mit sample_data
    matching_rows = diplotypes_df_filtered[
        diplotypes_df_filtered.apply(lambda row: all(row[key] == value for key, value in sample_data.items() if key in row), axis=1)
    ]

    results = match_and_sort(diplotypes_df, sample_data)
    # Ausgabe der übereinstimmenden Zeilen
    print(results)

    # # Evaluate matches
    # matches = evaluate_matches(diplotypes_df, sample_data, user_cnv)

    # # Print results
    # if matches:
    #     print_matches(diplotypes_df, matches)
    # else:
    #     print("No matches found.")


if __name__ == '__main__':
    main()
