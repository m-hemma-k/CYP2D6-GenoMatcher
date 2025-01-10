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

def evaluate_matches(diplotypes_df, data_input):
    """
    Match user-provided genetic input against a diplotypes DataFrame.

    Args:
        diplotypes_df (pd.DataFrame): DataFrame containing SNPs and genotypes.
        data_input (dict): Dictionary with SNP rsIDs as keys and allele values as strings.

    Returns:
        list: Genotype names from the DataFrame that match the user input.
    """
    # Create an empty list for matches
    matches = []

    # Convert data_input into a DataFrame
    user_df = pd.DataFrame(data_input, index=[0])

    for i, row in diplotypes_df.iterrows():
        try:
            # Check if all alleles match for the SNPs present in both user inputs and the DataFrame
            allele_match = all([
                set(user_df.loc[0, rsID].split('/')) == set(str(row[rsID]).split('/'))
                for rsID in data_input.keys()  # Check only keys in data_input
                if rsID in row.index  # Ensure rsID exists in the DataFrame
            ])

            # If there's a match, append the index or relevant column (e.g., 'Genotype') to the results
            if allele_match:
                matches.append(row['Genotype'])  # Add 'Genotype' to the match list

        except KeyError as e:
            # Handle cases where the data_input have an SNP not in the DataFrame
            print(f"KeyError: {e}. Check if all SNPs in data_input exist in the DataFrame.")

    return matches

def print_matches(diplotypes_df, matches):
    """
    Print matched genotypes sorted by ranking and include CNV values.

    Args:
        diplotypes_df (pd.DataFrame): DataFrame containing genotypes, rankings, and CNV values.
        matches (list): List of genotype names to filter and print.

    Returns:
        None
    """
    # Filter rows in diplotypes_df that match the Genotype in matches
    matches_df = diplotypes_df[diplotypes_df['Genotype'].isin(matches)]

    # Sort the filtered DataFrame by Ranking in descending order
    matches_sorted = matches_df.sort_values(by='Ranking', ascending=False)

    # Print the sorted matches with Ranking and CNV values
    print("Matches nach Wahrscheinlichkeit absteigend sortiert:")
    for _, row in matches_sorted.iterrows():
        genotype = row['Genotype']
        ranking = row['Ranking']
        cnv_value = row['CNV']
        print(f"100% Match: {genotype}, Ranking: {ranking}, CNV-Wert: {cnv_value}")

def main():
    """
    Main function to execute the genotype matching pipeline:
    1. Load the reference DataFrame containing diplotypes.
    2. Load and process VCF data to extract genotypes.
    3. Extract the CNV value from the VCF data.
    4. Compare user-provided genotypes against the reference DataFrame.
    5. Identify and rank matches based on likelihood and CNV values.
    6. Print the results sorted by ranking and CNV values.
    """

    # Step 1: Load the reference DataFrame containing diplotypes
    file_path = '../00-preprocessing_data/output/CYP2D6_2025-01-09.pkl'
    diplotypes_df = pd.read_pickle(file_path)

    # Step 2: Load the VCF file with sample genotype data
    sample_filepath = './input/Test_Input_1_1.vcf'
    vcf_data = read_vcf(sample_filepath)

    # Step 3: Extract the CNV value from the VCF data
    CNV = vcf_data.loc[vcf_data['ID'] == 'CYP2D6_CNV', 'PharmCAT']
    CNV_value = CNV.iloc[0] if not CNV.empty else None

    # Step 4: Extract genotypes for specific rsIDs from the VCF data
    sample_data = extract_rs_genotypes(diplotypes_df, vcf_data)

    # Step 5: Compare the extracted genotypes against the reference DataFrame
    results = evaluate_matches(diplotypes_df, sample_data)

    # Step 6: Print the matched genotypes sorted by ranking and CNV values
    print_matches(diplotypes_df, results)


if __name__ == '__main__':
    main()