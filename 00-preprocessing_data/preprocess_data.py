#!/usr/bin/env python3

import pandas as pd # type: ignore
from itertools import combinations_with_replacement
import re
import os
from datetime import datetime

def process_tsv_and_vcf(tsv_dir, vcf_dir, output_path):
    """
    Updates a TSV file with data from multiple VCF files.

    Args:
        tsv_path (str): Path to the input TSV file.
        vcf_dir (str): Directory containing VCF files.
        output_path (str): Path to save the updated TSV file.

    Returns:
        None: The updated TSV file is saved to the specified output path.
    """
    # Automatically find the TSV file in the directory
    tsv_file = next((file for file in os.listdir(tsv_dir) if file.endswith(".tsv")), None)
    if not tsv_file:
        raise FileNotFoundError("No TSV file found in the specified directory.")

    tsv_path = os.path.join(tsv_dir, tsv_file)

    # Load the TSV file and adjust header
    tsv_data = pd.read_csv(tsv_path, sep="\t", skiprows=1, header=0)

    # Ensure Variant Start and Variant Stop are integers
    tsv_data["Variant Start"] = pd.to_numeric(tsv_data["Variant Start"], errors="coerce").fillna(0).astype(int)
    tsv_data["Variant Stop"] = pd.to_numeric(tsv_data["Variant Stop"], errors="coerce").fillna(0).astype(int)

    # Filter rows containing 'insertion' or 'deletion'
    filtered_data = tsv_data[tsv_data["Type"].str.contains("insertion|deletion", case=False, na=False)]

    # Process each VCF file one by one
    for vcf_file in os.listdir(vcf_dir):
        if vcf_file.endswith(".vcf"):
            vcf_path = os.path.join(vcf_dir, vcf_file)

            with open(vcf_path, "r") as vcf:
                for line in vcf:
                    if line.startswith("#"):  # Skip header lines
                        continue
                    fields = line.strip().split("\t")
                    pos, ref, alt = int(fields[1]), fields[3], fields[4]

                    # Check if the POS matches any in the TSV file
                    for index, row in filtered_data.iterrows():
                        variant_start = row["Variant Start"]

                        # Adjust comparison for deletion
                        if "deletion" in row["Type"].lower():
                            match_pos = pos + 1
                        else:
                            match_pos = pos

                        if variant_start == match_pos:  # Match based on POS
                            # Update Variant Allele and Alternate Allele for insertion/deletion
                            tsv_data.at[index, "Reference Allele"] = ref
                            tsv_data.at[index, "Variant Allele"] = alt

    # Save the updated TSV data to a new file
    tsv_data.to_csv(output_path, sep="\t", index=False)

def filter_for_rsnumbers(pharmvardata, specific_rsIds):
    """
    Filter PharmVar data to include only the specified rsIDs and ensure haplotypes 
    'CYP2D6*5' and 'CYP2D6*1' are included at the beginning.

    Args:
        pharmvardata (DataFrame): Input DataFrame containing PharmVar data.
        specific_rsIds (list): List of specific rsIDs to filter.

    Returns:
        dict: Filtered haplotype data as a dictionary, with main haplotypes prioritized.
    """
    # Pre-filter the data to include only relevant rsIDs
    filtered_data = pharmvardata[pharmvardata['rsID'].isin(specific_rsIds)]

    # Convert Reference names and Variant names into dictionaries for fast lookups
    reference_dict = filtered_data.set_index('rsID')['Reference Allele'].to_dict()
    variant_dict = filtered_data.set_index(['Haplotype Name', 'rsID'])['Variant Allele'].to_dict()

    # Initialize the haplotype data structure
    haplotype_data = {}
    haplotypes = list(filtered_data['Haplotype Name'].unique())

    # Filter PharmVar for the specific rsIds and include *1 and *5
    # 'CYP2D5*1' would be missing, because it has no variant
    # 'CYP2D6*5' is not in the pharmvar tsv data, because it is a gene deletion. 
    # It will be available for analysis with CNV = 0 later on.
    # Filter PharmVar for the specific rsIds and include *1 and *5
    required_haplotypes = ['CYP2D6*5', 'CYP2D6*1']
    for haplotype in required_haplotypes:
        if haplotype not in haplotypes:
            haplotypes.insert(0, haplotype)  # Add missing haplotypes at the beginning

    # Iterate over each haplotype and rsID
    for haplotype_name in haplotypes:
        haplotype_data[haplotype_name] = {}
        for rsID in specific_rsIds:
            # Check for the variant value_name; if not found, use the reference value_name
            variant_name = variant_dict.get((haplotype_name, rsID))
            reference_name = reference_dict.get(rsID, '-')
            haplotype_data[haplotype_name][rsID] = variant_name if variant_name is not None else reference_name

    # check if a subvariant is equal to the main variant, if yes, remove the subvariant
    for haplotype_name in list(haplotype_data.keys()):
    # Check if the haplotype value_name has a subvariant (contains a '.')
        if '.' in haplotype_name:
            # Extract the main variant value_name (everything before the first '.')
            main_variant = haplotype_name.split('.')[0]

            # Check if the main variant exists in haplotype_data
            if main_variant in haplotype_data:
                # Compare the key-value pairs of the main variant and the subvariant
                if haplotype_data[haplotype_name] == haplotype_data[main_variant]:
                    # If identical, remove the subvariant
                    del haplotype_data[haplotype_name]

    return haplotype_data

def special_combinations(haplotype_data, special_combinations, specific_rsIds):
    """
    Generate combined haplotypes for special combinations and update the haplotype data.

    Args:
        haplotype_data (dict): Dictionary containing haplotype data.
        special_combinations (set): Set of tuples defining special haplotype combinations.
        specific_rsIds (list): List of specific rsIDs to process.

    Returns:
        dict: Updated haplotype data including special combinations.
    """
    # Iterate over each special combination
    for combination in special_combinations:
        combined_name = "+".join(combination)  # Create a combined haplotype value_name
        haplotype_data[combined_name] = {}    # Add new combined haplotype

        # Process each rsID in the specific list
        for rsID in specific_rsIds:
            # Collect names for this rsID from the parts of the combination
            names = [
                haplotype_data[part][rsID]
                for part in combination
                if part in haplotype_data and rsID in haplotype_data[part]
            ]

            # Determine the resulting value_name for the combination
            if not names:  # No names found
                name_value = '-'
            elif len(set(names)) == 1:  # All names are the same
                name_value = names[0]
            else:  # Combine differing names with "/"
                name_value = '/'.join(sorted(set(names)))

            # Update haplotype_data directly
            haplotype_data[combined_name][rsID] = name_value

    return haplotype_data

def add_cnv_values_ex9(haplotype_data, cnv_exceptions, special_combinations_cnv):
    """
    Add CNV values analysed in exon 9 to haplotypes, taking into account exceptions and special combinations.

    Args:
        haplotype_data (dict): Dictionary containing haplotype data.
        cnv_exceptions (list): List of haplotypes with CNV set to 0.
        special_combinations_cnv (dict): Dictionary of special CNV combinations with predefined values.

    Returns:
        dict: Updated haplotype data with CNV values assigned.
    """
    # Add the CNV column to the haplotype_data dictionary
    for haplotype_name in haplotype_data:
        # Default CNV value is set to 1 for each haplotype
        cnv_value = 1

        # Check if the haplotype is in the exceptions list and set CNV to 0 if it is
        if haplotype_name in cnv_exceptions:
            cnv_value = 0

        # Assign the CNV value to the current haplotype
        haplotype_data[haplotype_name]["CNV"] = cnv_value

        # If the haplotype is part of a special combination, use the predefined CNV value
        if haplotype_name in special_combinations_cnv:
            # Use the predefined CNV value for the special combination
            cnv_value = special_combinations_cnv[haplotype_name]

            # Set the CNV value for the special combination
            haplotype_data[haplotype_name]["CNV"] = cnv_value

    return haplotype_data

def cnv_variations(haplotype_data, cnv_duplications, cnv_triplication):
    """
    Create CNV duplications and triplications for specified haplotypes.

    Args:
        haplotype_data (dict): Dictionary containing haplotype data.
        CNV_DEFINITIONS_2 (list): List of haplotypes of gene duplication to get CNV = 2.
        CNV_DEFINITIONS_3 (list): List of haplotypes of gene triplication to get CNV = 3.

    Returns:
        dict: Updated haplotype data including CNV duplications and triplications.
    """
    # Duplicate haplotypes for CNV = 2
    for haplotype_name in cnv_duplications:
        if haplotype_name in haplotype_data:
            duplicated_entry = haplotype_data[haplotype_name].copy()
            duplicated_entry["CNV"] = 2
            haplotype_data[f"{haplotype_name}x2"] = duplicated_entry

    # Duplicate haplotypes for CNV = 3
    for haplotype_name in cnv_triplication:
        if haplotype_name in haplotype_data:
            duplicated_entry = haplotype_data[haplotype_name].copy()
            duplicated_entry["CNV"] = 3
            haplotype_data[f"{haplotype_name}x3"] = duplicated_entry

    return haplotype_data

def add_ranking(haplotype_data, ranking):
    """
    Assign ranking values (Value: 0, 1, or 2) to haplotypes based on predefined criteria.

    Args:
        haplotype_data (dict): Dictionary containing haplotype data.
        ranking (dict): Dictionary defining 'Top Tier' and '2nd Tier' haplotypes.

    Returns:
        dict: Updated haplotype data with ranking values assigned.
    """
    # Create a dictionary to store the Tier values
    tier_dict = {}

    for haplotype_name in haplotype_data.keys():
        if haplotype_name in ranking.get('Top Tier', []):
            tier_dict[haplotype_name] = 2
        elif haplotype_name.endswith("x2"):
            tier_dict[haplotype_name] = 2
        elif haplotype_name.endswith("x3"):
            tier_dict[haplotype_name] = 2
        elif haplotype_name in ranking.get('2nd Tier', []):
            tier_dict[haplotype_name] = 1
        elif '+' in haplotype_name:
            tier_dict[haplotype_name] = 1
        else:
            tier_dict[haplotype_name] = 0

    # Update the haplotype_data dictionary with the assigned Tier values
    for haplotype_name in haplotype_data:
        haplotype_data[haplotype_name]["Ranking"] = tier_dict.get(haplotype_name, 0)

    return haplotype_data

def extract_numeric_value(haplotype_name):
    """
    Extract the numeric value following 'CYP2D6*' in a haplotype name.

    Args:
        haplotype_name (str): Haplotype name, e.g., 'CYP2D6*41x3'.

    Returns:
        int: Extracted numeric value. Returns 0 if no numeric value is found.
    """
    match = re.search(r'\*([0-9]+)', haplotype_name)
    return int(match.group(1)) if match else 0

def pair_haplotypes(haplotype_data):
    """
    Generate all unique pairs of haplotypes (including self-pairing), combining their 
    CNV, Ranking, and rsID values.

    Args:
        haplotype_data (dict): Dictionary containing haplotype data.

    Returns:
        dict: New dictionary with paired haplotypes and their combined values.
    """
    # List of all unique haplotypes
    haplotypes_list = list(haplotype_data.keys())

    # Dictionary to store the combined results
    combinations_dict = {}

    # Iterate over all combinations (including self-pairing)
    for haplotype1, haplotype2 in combinations_with_replacement(haplotypes_list, 2):
        # Sort numerically based on the number after 'CYP2D6'
        sorted_haplotypes = sorted([haplotype1, haplotype2], key=extract_numeric_value)
        combined_key = f"{sorted_haplotypes[0]}/{sorted_haplotypes[1]}"

        if combined_key not in combinations_dict:
            combinations_dict[combined_key] = {}

        # Process each attribute (rs values, CNV, Ranking)
        for value_name in set(haplotype_data[haplotype1].keys()).union(haplotype_data[haplotype2].keys()):
            value1 = haplotype_data[haplotype1].get(value_name)
            value2 = haplotype_data[haplotype2].get(value_name)

            # Special handling for CNV and Ranking
            if value_name in ["CNV", "Ranking"]:
                combined_value = value1 + value2
            else:
                combined_value = value1 if value1 == value2 else f"{value1}/{value2}"

            # Store the combined value
            combinations_dict[combined_key][value_name] = combined_value

    return combinations_dict

def combinations_dict_to_dataframe(combinations_dict):
    """
    Convert a dictionary of haplotype combinations into a Pandas DataFrame.

    Args:
        combinations_dict (dict): Dictionary of haplotype combinations with combined values.

    Returns:
        DataFrame: Pandas DataFrame with haplotype combinations and associated values.
    """
    # Convert dictionary to DataFrame
    dataframe = pd.DataFrame.from_dict(combinations_dict, orient='index')

    # Reset index to make combined keys a column
    dataframe.reset_index(inplace=True)
    dataframe.rename(columns={'index': 'Genotype'}, inplace=True)

    return dataframe

def save_dataframe_to_pickle(dataframe, save_directory):
    """
    Save a Pandas DataFrame to a pickle file in the specified directory.

    Args:
        dataframe (DataFrame): The DataFrame to save.
        save_directory (str): Path to the directory where the file will be saved.

    Returns:
        str: Path to the saved pickle file.
    """
    # Ensure the save directory exists
    os.makedirs(save_directory, exist_ok=True)

    # Get the current date and time for the file value_name
    current_date = datetime.now().strftime("%Y-%m-%d")
    file_name = f"CYP2D6_{current_date}.pkl"
    file_path = os.path.join(save_directory, file_name)

    # Save the DataFrame to the specified file path
    dataframe.to_pickle(file_path)
    print(f"File saved: {file_path}")

    return file_path

# Hauptprogramm
def main():
    """
    Main function to execute the processing pipeline:
    1. Rewrite Reference Allele and Stop with ALT and REF from vcf files.
    2. Load PharmVar data and specific rsIDs.
    3. Filter and preprocess haplotype data.
    4. Generate special combinations.
    5. Assign CNV Exon 9 values, duplications, triplications and rankings.
    6. Pair haplotypes and convert to a DataFrame.
    7. Save the final DataFrame to a pickle file.
    """

    # Rewrite data with vcf files
    tsv_path = './input/RefSeqGene'
    vcf_dir = './input/RefSeqGene'
    output_path = './input/CYP2D6.tsv'
    process_tsv_and_vcf(tsv_path, vcf_dir, output_path)

    # Read in data from PharmVar
    file_path = output_path
    pharm_var_data = pd.read_csv(file_path, sep='\t', header=0, comment='#')

    # Read in specified rsIDs
    rsnumbers_data = './input/rs_numbers.txt'
    with open(rsnumbers_data, 'r') as file:
        rs_ids = [line.strip() for line in file if line.strip()]

    # Filter PharmVar for the specific rsIds and include *1 and *5
    # 'CYP2D5*1' would be missing, because it has no variant
    # 'CYP2D6*5' is not in the pharmvar tsv data, because it is a gene deletion. 
    # It will be available for analysis with CNV = 0 later on.
    # Remove Subvariants if they are equal to their main variant.
    pharmvar_filtered_data = filter_for_rsnumbers(pharm_var_data, rs_ids)

    # add hybrid genes
    hybrid_genes = {
        ("CYP2D6*1", "CYP2D6*38"),
        ("CYP2D6*4.013", "CYP2D6*4"),
        ("CYP2D6*13", "CYP2D6*1"),
        ("CYP2D6*13", "CYP2D6*2"),
        ("CYP2D6*13", "CYP2D6*68", "CYP2D6*4"),
        ("CYP2D6*17", "CYP2D6*17"),
        ("CYP2D6*36", "CYP2D6*10"),
        ("CYP2D6*36", "CYP2D6*10.007"),
        ("CYP2D6*36.004", "CYP2D6*10.002"),
        ("CYP2D6*57", "CYP2D6*10"),
        ("CYP2D6*68", "CYP2D6*2"),
        ("CYP2D6*68", "CYP2D6*4"),
        ("CYP2D6*1", "CYP2D6*90")
    }
    all_haplotypes_data = special_combinations(pharmvar_filtered_data, hybrid_genes, rs_ids)
    
    # Assign CNV values to haplotypes
    cnv_exon_9_zero = [
        'CYP2D6*4.013', 'CYP2D6*4.031', 'CYP2D6*36', 'CYP2D6*36.001', 
        'CYP2D6*36.002', 'CYP2D6*36.003', 'CYP2D6*36.004', 'CYP2D6*36.005',
        'CYP2D6*83', 'CYP2D6*83.001', 'CYP2D6*83.0012', 'CYP2D6*83.0013',
        'CYP2D6*141', 'CYP2D6*141.001', 'CYP2D6*5'
    ]
    cnv_exon9_hybrid_genes = {
        "CYP2D6*1+CYP2D6*38": 2,
        "CYP2D6*4.013+CYP2D6*4": 1,
        "CYP2D6*13+CYP2D6*1": 2,
        "CYP2D6*13+CYP2D6*2": 2,
        "CYP2D6*13+CYP2D6*68+CYP2D6*4": 3,
        "CYP2D6*17+CYP2D6*17": 2,
        "CYP2D6*36+CYP2D6*10": 1,
        "CYP2D6*36+CYP2D6*10.007": 1,
        "CYP2D6*36.004+CYP2D6*10.002": 1,
        "CYP2D6*57+CYP2D6*10": 2,
        "CYP2D6*68+CYP2D6*2": 2,
        "CYP2D6*68+CYP2D6*4": 2,
        "CYP2D6*1+CYP2D6*90": 2,
        "CYP2D6*28.001+CYP2D6*28.003": 2
    }
    haplotypes_data_cnv_ex9 = add_cnv_values_ex9(all_haplotypes_data, cnv_exon_9_zero, cnv_exon9_hybrid_genes)

    # include duplications and triplications
    cnv_duplications = [
        'CYP2D6*1', 'CYP2D6*2', 'CYP2D6*3', 'CYP2D6*4', 'CYP2D6*4.013', 'CYP2D6*6',
        'CYP2D6*9', 'CYP2D6*10', 'CYP2D6*17', 'CYP2D6*27.002', 'CYP2D6*29',
        'CYP2D6*35', 'CYP2D6*41', 'CYP2D6*43', 'CYP2D6*45', 'CYP2D6*146.001'
    ]
    cnv_triplication = ['CYP2D6*1', 'CYP2D6*2', 'CYP2D6*4', 'CYP2D6*41']
    data_with_cnv = cnv_variations(haplotypes_data_cnv_ex9, cnv_duplications, cnv_triplication)

    # Add a ranking according to a Nature Paper - CYP2D6 Overview: Alelle and Phenotype Frequencies, Megan Kane, PhD., October 15, 2021.
    ranking = {
        'Top Tier': ['CYP2D6*1', 'CYP2D6*2', 'CYP2D6*3', 'CYP2D6*4', 'CYP2D6*5', 'CYP2D6*6', 'CYP2D6*9', 'CYP2D6*10', 'CYP2D6*17', 'CYP2D6*29', 'CYP2D6*41'],
        '2nd Tier': ['CYP2D6*7', 'CYP2D6*8', 'CYP2D6*12', 'CYP2D6*14', 'CYP2D6*15', 'CYP2D6*21', 'CYP2D6*31', 'CYP2D6*40', 'CYP2D6*42', 'CYP2D6*49', 'CYP2D6*56', 'CYP2D6*59']
    } # Multiplications are also 'Top Tier', hybrid genes are '2nd Tier'.
    haplotypes_dic = add_ranking(data_with_cnv, ranking)

    # Pair each haplotype
    combinations = pair_haplotypes(haplotypes_dic)

    # Convert dict to pandaframe
    pf_combinations = combinations_dict_to_dataframe(combinations)
    # print(pf_combinations)

    # Save dataframe to pkl
    filepath = './output/'
    save_dataframe_to_pickle(pf_combinations, filepath)

if __name__ == '__main__':
    main()