
---

# CYP2D6 GenoMatcher
#### *For Research Use Only*
This repository contains Python scripts designed to process haplotype data from PharmVar and analyze input variant data formatted as a VCF (Variant Call Format) file. The tool focuses specifically on the **CYP2D6 gene**, providing an efficient workflow for filtering, combining, and analyzing haplotypes to ultimately determine genotypes.

The current version of the program is designed to extract the Copy Number Variation (CNV) value specifically for exon 9. This functionality is based on the Thermo Fisher Scientific assay Hs00010001_cn. To ensure proper usage and interpretation of the assay data, users are advised to carefully read the instructions provided on the official Thermo Fisher website: [Thermo Fisher Assay Hs00010001_cn](https://www.thermofisher.com/order/genome-database/details/copy-number/Hs00010001_cn).

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

   Assigns rankings to haplotypes based on scientific categorizations, as outlined by Megan Kane, PhD ([NCBI Bookshelf, 2021](https://www.ncbi.nlm.nih.gov/books/NBK574601/):

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

3. **Extracting CNV values**

   - Extract CNV values from a txt file.

3. **User Input Matching**

   - Compares user-provided genetic input against a reference DataFrame of haplotypes.

4. **Match Evaluation**

   - Matches user input with haplotypes, sorts them by ranking, and includes CNV values.

5. **Output**

   - Prints the genotype match directly to a text file (e.g., `CYP2D6*1/CYP2D6*1`).

---

### Installation

---

## **Prepare Input Files**

Before running the installation script, ensure the required input files are in their predefined locations. 
These files are already pre-saved but can be modified if needed by following these steps:

1. **PharmVar Data**:

   - Navigate to [PharmVar CYP2D6](https://www.pharmvar.org/gene/CYP2D6).
   - Download the haplotype data (`Download Gene Data`).
   - Extract the contents and place the `CYP2D6-6.2` folder in `00-preprocessing_data/input`.

2. **rs\_numbers.txt**:

   - Ensure `rs_numbers.txt` is located in the `00-preprocessing_data/input` directory.
   - You can customize the predefined list as needed.

   A list of `rsIDs`, e.g.:

   ```
   rs35742686
   rs3892097
   ```

---

## **Step 2: Run the Installation Script**

The installation script is located in the `00-preprocessing_data` directory.

Execute the `install.sh` script in the terminal:

```bash
./install.sh
```

The `install.sh` script performs the following actions:

- **Dependency Installation**: Ensures that the required Python dependencies (e.g., `pandas`) are installed.
- **Preprocessing Execution**: Runs the `preprocess_data.py` script to prepare initial data.
- **Genotype Script Setup**: Makes the `genotype_CYP2D6.py` script executable.
- **User-Friendly Shortcut Creation**: Creates a `.bat` file (`Run_Genotype.bat`) in the `01-genotype_CYP2D6` folder for running the genotype evaluation script without requiring terminal interaction.

### **Step 3: Verify Setup**

After running `install.sh`, the environment should be ready for analysis. Proceed to the next section to run the genotype evaluation.

---

## **Usage**

### **1. Preprocessing Script**

1. **Sample Files**:

   For the first usage two files are already presaved in the `01-genotype_CYP2D6/input` folder:
   - 
   - 

   For further usage:
   - Place your sample `.vcf` file in the `01-genotype_CYP2D6/input` folder.
   - Place your sample `.txt` file in the `01-genotype_CYP2D6/input` folder.
   - The `.vcf` and `.txt` file MUST have the same name, otherwise they cannot be processed.

   Here is an example of a `.vcf` file:

   ```
   ##fileformat=VCFv4.2
   ##source=PharmCAT allele definitions
   ##fileDate=2023-02-03T08:20:52.044-08:00
   ##contig=<ID=chr1,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr2,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr4,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr6,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr7,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr10,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr12,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr13,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr16,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr19,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chr22,assembly=GRCh38.p13,species="Homo sapiens">
   ##contig=<ID=chrX,assembly=GRCh38.p13,species="Homo sapiens">
   ##FILTER=<ID=PASS,Description="All filters passed">
   ##INFO=<ID=PX,Number=.,Type=String,Description="Gene">
   ##INFO=<ID=POI,Number=0,Type=Flag,Description="Position of Interest but not part of an allele definition">
   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	PharmCAT
   chr22	42130692	rs1065852	G	A	.	PASS	PX=CYP2D6	GT	0/0
   chr22	42132375	rs1080985	G	C	.	PASS	PX=CYP2D6	GT	0/0
   ```

   The `.txt` file must adhere to this format: the first and second words should be separated by a tab:

   ```
   Exon9 2
   ```

### **2. Genotype Evaluation Script**

To evaluate genotypes, use the shortcut created during the installation:

1. Navigate to the `01-genotype_CYP2D6` folder.
2. Double-click the `Run_Genotype.bat` file to execute the script.

Alternatively, you can run the script directly from the terminal:

```bash
python 01-genotype_CYP2D6/genotype_CYP2D6.py
```

### **3. Output**

The output will be saved as a `.txt` file in the `01-genotype_CYP2D6/output` folder, containing the best genotype result (e.g., `CYP2D6*1/CYP2D6*1`).
The input sample files will be moved to the folder `01-genotype_CYP2D6/input/processed_data`.

---

## **Notes**

### **Error Handling**

1. If the `install.sh` script reports missing files or directories:

   - Verify that all input files are placed in their respective locations.
   - Check that required folder `00/-preprocessong_data/input/CYP2D6-6.2/GRCh38` exists.

2. If the `.bat` file does not execute:

   - Ensure that Python is installed and added to the system PATH.
   - Confirm that the correct Python version is being used (3.12 or later).

### **System Requirements**

- Python 3.12 or later
- Required Python packages (installed automatically by `install.sh`):
  - `pandas >= 2.1.2`

---

## **License**

This repository is licensed under the **MIT License**. Please refer to the `LICENSE` file for details.

---

## **Citation**

If you use this script or any part of it for **publication or presentation**, please cite:

**Hemma Kargl, 2024** *"CYP2D6 GenoMatcher."*

### References:

1. **CYP2D6 Overview: Genotype and Phenotype Frequencies**   Source: Megan Kane, PhD, October 15, 2021   Available at: [NCBI Bookshelf - CYP2D6](https://www.ncbi.nlm.nih.gov/books/NBK574601/)

2. **PharmVar Consortium**   PharmVar Gene: CYP2D6   Available at: [https://www.pharmvar.org/gene/CYP2D6](https://www.pharmvar.org/gene/CYP2D6)   Accessed: 11/01/2025

3. Gaedigk A, Ingelman-Sundberg M, Miller NA, et al.   The Pharmacogene Variation (PharmVar) Consortium: Incorporation of the Human Cytochrome P450 (CYP) Allele Nomenclature Database.   *Clin Pharmacol Ther.* 2018;103(3):399-401.   DOI: [10.1002/cpt.910](https://doi.org/10.1002/cpt.910)

4. Thermo Fisher Scientific - Genome Database
Copy Number Assay: Hs00010001_cn
Available at: https://www.thermofisher.com/order/genome-database/details/copy-number/Hs00010001_cn
Accessed: 11/01/2025

---

## **Contact**

For any questions, issues, or suggestions, please contact here or on LinkedIn [[www.linkedin.com/in/hemma-kargl](http://www.linkedin.com/in/hemma-kargl)].
