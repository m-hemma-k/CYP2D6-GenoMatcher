---

# CYP2D6 Haplotype Processing and Analysis

This repository contains a Python script designed to preprocess and analyze haplotype data, specifically focusing on the **CYP2D6 gene**. The script provides an efficient workflow for filtering, combining, and analyzing haplotypes, enabling researchers to study copy number variations (CNV), hybrid genes, and rankings of CYP2D6 haplotypes. The final processed data is output as a **Pandas DataFrame** and saved in a `.pkl` format for further analysis.

---

## Features

1. **Rewrite Reference Allele and Stop**
   - Updates Variant Allele and Alternate Allele in the TSV file using ALT and REF values from the VCF files.

2. **Filtering and Preprocessing**  
   - Filters PharmVar data to include only relevant `rsIDs`.
   - Ensures that haplotypes like `CYP2D6*1` and `CYP2D6*5` are included in the analysis.

3. **Handling Special Combinations**  
   - Processes hybrid haplotypes (e.g., `CYP2D6*1 + CYP2D6*38`).
   - Combines alleles for specified `rsIDs` into a unified dataset.

4. **Copy Number Variation (CNV)**  
   - Generates duplications (`x2`) and triplications (`x3`) for specified haplotypes.
   - Assigns predefined CNV values for specific hybrid haplotypes.

5. **Haplotype Ranking**  
   - Assigns rankings to haplotypes based on scientific categorizations:
     - **Top Tier (2)**: Common and impactful haplotypes.
     - **Second Tier (1)**: Less impactful haplotypes.
     - **Others (0)**: Remaining haplotypes.

6. **Unique Pair Generation**  
   - Creates unique pairings of haplotypes, including self-pairings.
   - Combines CNV values, rankings, and allele data for each pairing.

7. **Output**  
   - Saves the processed data as a `.pkl` file in a user-specified directory.
   - Outputs are ready for direct use in downstream analysis or visualization.

---

## Installation

1. Clone this repository:
   ```bash
   git clone git@github.com:m-hemma-k/CYP2D6.vcf.git
   ```

2. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Prepare your input files:
   ```
   The necessary files and folders are saved in the `input` directory by default but can be modified if needed.
   - **PharmVar Data**: Contains a haplotype data folder downloaded from PharmVar.
                        - Navigate to the PharmVar CYP2D6 page: [PharmVar CYP2D6](https://www.pharmvar.org/gene/CYP2D6).
                        - Click on `Download Gene Data` to download the dataset.
                        - Extract the downloaded file and move the `RefSeqGene` folder into the `input` directory.
   - **rs_numbers.txt**: A text file listing specific `rsIDs` to include in the analysis.
   ```

---

## Usage

1. Run the script:
   ```bash
   python script_name.py
   ```

2. The output will be saved in the specified output directory as a `.pkl` file.

---

## Input and Output

### Input
1. **PharmVar Data**:

   **1.1 `.tsv` File**
   The TSV file must follow this structure:
   ```
   #version=pharmvar-6.1.8
   Haplotype Name	Gene	rsID	ReferenceSequence	Variant Start	Variant Stop	Reference Allele	Variant Allele	Type
   CYP2D6*1	CYP2D6	REFERENCE	.	.	.	.	.	.
   CYP2D6*1.001	CYP2D6	REFERENCE	.	.	.	.	.	.
   CYP2D6*1.002	CYP2D6	rs28371732	NG_008376.4	8848	8848	G	A	substitution
   CYP2D6*1.003	CYP2D6	rs150163869	NG_008376.4	6998	6998	C	T	substitution
   CYP2D6*1.004	CYP2D6	rs28371718	NG_008376.4	7595	7595	C	A	substitution
   CYP2D6*1.005	CYP2D6	rs111606937	NG_008376.4	6889	6889	T	C	substitution
   CYP2D6*1.006	CYP2D6	rs28371699	NG_008376.4	5329	5329	G	T	substitution
   ```

   **1.2 `.vcf` Files**
   The VCF files should contain haplotype data in the following format:
   ```
   ##fileformat=VCFv4.1
   ##fileDate=20241216
   ##contig=<ID=NG_008376.4,length=11312>
   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   ##INFO=<ID=VI,Number=1,Type=String,Description="Variant impact">
   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
   NG_008376.4	4976	rs28371695	G	GG	.	.	.
   NG_008376.4	7780	rs28371723	T	A	.	.	.
   NG_008376.4	9413	rs202032066	CAG	C	.	.	.
   ```

3. **rs_numbers.txt**: A list of `rsIDs` to filter, e.g.:
   ```
   rs35742686
   rs3892097
   ```

### Output
The processed data is saved as a `.pkl` file. Example structure of the resulting DataFrame:
```
   Genotype                 rs5030655   rs28371725   ...   Ranking   CNV
   CYP2D6*1/CYP2D6*1        T           G            ...   2         2
   CYP2D6*1/CYP2D6*5        T           G            ...   2         2
   CYP2D6*1/CYP2D6*1.032    T           G            ...   2         2
   CYP2D6*1/CYP2D6*2        T           G            ...   2         2
   CYP2D6*1/CYP2D6*2.001    T           G            ...   2         2
   ...                      ...         ...          ...   ...       ...
   CYP2D6*2x3/CYP2D6*4x3    T           G            ...   2         5
   CYP2D6*2x3/CYP2D6*41x3   T           G/A          ...   2         5
   CYP2D6*4x3/CYP2D6*4x3    T           G            ...   2         6
   CYP2D6*4x3/CYP2D6*41x3   T           G/A          ...   2         6
   CYP2D6*41x3/CYP2D6*41x3  T           A            ...   2         6

```

---

## License

This repository is licensed under the **MIT License**.  
Please refer to the `LICENSE` file for details.

---

## Citation
If you use this script or any part of it for **publication or presentation**, please cite:

**Hemma Kargl, 2024**  
*"CYP2D6 Haplotype Processing and Analysis Script."*

Additionally, the haplotype data used in this script originates from the **PharmVar Consortium**. Please acknowledge and cite the PharmVar database as follows:

**PharmVar Consortium**  
PharmVar Gene: CYP2D6  
Available at: [https://www.pharmvar.org/gene/CYP2D6](https://www.pharmvar.org/gene/CYP2D6)  
Accessed: [insert date of access]
 
Gaedigk A, Ingelman-Sundberg M, Miller NA, et al. The Pharmacogene Variation (PharmVar) Consortium: Incorporation of the Human Cytochrome P450 (CYP) Allele Nomenclature Database. *Clin Pharmacol Ther.* 2018;103(3):399-401. doi:[10.1002/cpt.910](https://doi.org/10.1002/cpt.910)

**CYP2D6 Overview: Genotype and Phenotype Frequencies**
Source: Megan Kane, PhD., October 15, 2021.
Available at: https://www.ncbi.nlm.nih.gov/books/NBK574601/

---

## Contact

For any questions, issues, or suggestions, please contact here or on LinkedIn [www.linkedin.com/in/hemma-kargl].