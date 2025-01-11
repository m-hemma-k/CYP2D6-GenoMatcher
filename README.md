
---

# CYP2D6 GenoMatcher
***For Research Use Only***
This repository contains Python scripts designed to process haplotype data from PharmVar and analyze input variant data formatted as a VCF (Variant Call Format) file. The tool focuses specifically on the **CYP2D6 gene**, providing an efficient workflow for filtering, combining, and analyzing haplotypes to ultimately determine genotypes.

---

## Features

### Preprocessing Script:

This script should only be run when updating your desired `rsIDs` or when new PharmVar data becomes available.

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

   Assigns rankings to haplotypes based on scientific categorizations, as outlined by Megan Kane, PhD ([NCBI Bookshelf, 2021](https://www.ncbi.nlm.nih.gov/books/NBK574601/)):

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

### Main Script (Genotype Evaluation Script):

This is the main script to run for evaluating genotypes.

1. **Input VCF Files**

   - Reads VCF files and processes them into a DataFrame while ignoring metadata lines.

2. **Extracting rsID Genotypes**

   - Extracts rsIDs from the VCF file and matches them with a reference DataFrame.

3. **User Input Matching**

   - Compares user-provided genetic input against a reference DataFrame of haplotypes.

4. **Match Evaluation**

   - Matches user input with haplotypes, sorts them by ranking, and includes CNV values.

5. **Output**

   - Prints the genotype match directly to a text file (e.g., `CYP2D6*1/CYP2D6*1`).

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
   Place the required input files into the appropriate folders:
   - **PharmVar Data**: Contains a haplotype data folder downloaded from PharmVar.
                        - Navigate to the PharmVar CYP2D6 page: [PharmVar CYP2D6](https://www.pharmvar.org/gene/CYP2D6).
                        - Click on `Download Gene Data` to download the dataset.
                        - Extract the downloaded file and move the `RefSeqGene` folder into the `input` directory.
   - **VCF Files**: Place your sample VCF files into the `input` folder within the second folder.
   - **rs_numbers.txt**: A text file listing specific `rsIDs` to include in the analysis.
   ```

---

## Usage

1. If needed, run the preprocessing script to update `rsIDs` or PharmVar data:

   ```bash
   python preprocessing_script.py
   ```

2. To evaluate genotypes, run the primary script:

   ```bash
   python genotype_CYP2D6.py
   ```

3. The output will be saved as a `.txt` file with the best genotype result (e.g., `CYP2D6*1/CYP2D6*1`).

---

## Input and Output

### Input

1. **PharmVar Data**:

   **1.1 `.tsv` File**

   ```
   #version=pharmvar-6.1.8
   Haplotype Name	Gene	rsID	ReferenceSequence	Variant Start	Variant Stop	Reference Allele	Variant Allele	Type
   CYP2D6*1	CYP2D6	REFERENCE	.	.	.	.	.	.
   CYP2D6*1.001	CYP2D6	REFERENCE	.	.	.	.	.	.
   CYP2D6*1.002	CYP2D6	rs28371732	NG_008376.4	8848	8848	G	A	substitution
   ```

   **1.2 `.vcf` Files**

   ```
   ##fileformat=VCFv4.1
   ##fileDate=20241216
   ##contig=<ID=NG_008376.4,length=11312>
   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   ##INFO=<ID=VI,Number=1,Type=String,Description="Variant impact">
   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
   NG_008376.4	4976	rs28371695	G	GG	.	.	.
   ```

2. **rs_numbers.txt**: A list of `rsIDs` to filter, e.g.:

   ```
   rs35742686
   rs3892097
   ```

### Output

The processed data is saved as a `.txt` file containing only the best genotype result (e.g., `CYP2D6*1/CYP2D6*1`).

---

## License

This repository is licensed under the **MIT License**.Please refer to the `LICENSE` file for details.

---

## Citation

If you use this script or any part of it for **publication or presentation**, please cite:

**Hemma Kargl, 2024** *"CYP2D6 GenoMatcher."*

### References:

1. **CYP2D6 Overview: Genotype and Phenotype Frequencies**   Source: Megan Kane, PhD, October 15, 2021   Available at: [NCBI Bookshelf - CYP2D6](https://www.ncbi.nlm.nih.gov/books/NBK574601/)

2. **PharmVar Consortium**   PharmVar Gene: CYP2D6   Available at: [https://www.pharmvar.org/gene/CYP2D6](https://www.pharmvar.org/gene/CYP2D6)   Accessed: 11/01/2025

3. Gaedigk A, Ingelman-Sundberg M, Miller NA, et al.   The Pharmacogene Variation (PharmVar) Consortium: Incorporation of the Human Cytochrome P450 (CYP) Allele Nomenclature Database.   *Clin Pharmacol Ther.* 2018;103(3):399-401.   DOI: [10.1002/cpt.910](https://doi.org/10.1002/cpt.910)

---

## Contact

For any questions, issues, or suggestions, please contact here or on LinkedIn [[www.linkedin.com/in/hemma-kargl](http://www.linkedin.com/in/hemma-kargl)].
