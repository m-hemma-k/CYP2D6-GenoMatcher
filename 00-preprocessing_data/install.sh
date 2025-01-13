#!/bin/bash

# Installation and Initialization Script for CYP2D6 GenoMatcher

# Define paths
SUBFOLDER="01-genotype_CYP2D6"
PREPROCESSING_SCRIPT="preprocess_data.py"
GENOTYPE_SCRIPT="../$SUBFOLDER/genotype_CYP2D6.py"
DESKTOP_FILE="Run_Genotype.desktop"

# Use Python explicitly
PYTHON_EXEC="python"  # Use 'python' since it's associated with Python 3.12

# Step 1: Check for pandas
echo "Checking if pandas is installed..."
if ! "$PYTHON_EXEC" -c "import pandas" &>/dev/null; then
    echo "pandas is not installed. Installing pandas..."
    "$PYTHON_EXEC" -m pip install pandas>=2.1.2 || { echo "Failed to install pandas."; exit 1; }
else
    echo "pandas is already installed."
fi

# Step 2: Run the preprocessing script
if [ -f "$PREPROCESSING_SCRIPT" ]; then
    echo "Running the preprocessing script to initialize data..."
    "$PYTHON_EXEC" "$PREPROCESSING_SCRIPT" || { echo "Failed to run preprocess_data.py."; exit 1; }
else
    echo "Error: preprocess_data.py not found."
    exit 1
fi

# Step 3: Make genotype.py executable
if [ -f "$GENOTYPE_SCRIPT" ]; then
    echo "Making $GENOTYPE_SCRIPT executable..."
    chmod +x "$GENOTYPE_SCRIPT" || { echo "Failed to make $GENOTYPE_SCRIPT executable."; exit 1; }
else
    echo "Error: genotype_CYP2D6.py not found in 01-genotype_CYP2D6."
    exit 1
fi

# Step 4: Create a .bat file to run genotype.py on Windows
BAT_FILE="../01-genotype_CYP2D6/Run_Genotype.bat"

echo "Creating a .bat file for $GENOTYPE_SCRIPT in the subfolder..."
cat <<EOF > $BAT_FILE
@echo off
echo Running CYP2D6 Genotype Script...
cd /d "$(dirname "$(pwd)/$GENOTYPE_SCRIPT")"
python "$(basename "$GENOTYPE_SCRIPT")"
pause
EOF

# Set permissions for the .bat file (optional on Windows)
chmod +x "$BAT_FILE"

echo "Installation completed successfully!"
echo "You can now double-click the 'Run_Genotype' file to run genotype.py."