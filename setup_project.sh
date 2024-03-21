#!/bin/bash

# Define the desired project name and environment name
DESIRED_PROJECT_NAME="aDNA_Comparative_Analysis"
ENV_NAME="aDNA_env"

# Get the name of the current directory
CURRENT_DIR_NAME=$(basename "$(pwd)")

# Check if the current directory is the desired project name
if [ "$CURRENT_DIR_NAME" != "$DESIRED_PROJECT_NAME" ]; then
    echo "Current directory is not '$DESIRED_PROJECT_NAME'. Creating the project directory..."
    mkdir -p "$DESIRED_PROJECT_NAME"
    cd "$DESIRED_PROJECT_NAME"
fi

# Set up directory structure with checks for existing directories
echo "Setting up project structure for $DESIRED_PROJECT_NAME..."

# Function to create directory if it doesn't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        mkdir -p "$1"
        echo "Created directory $1"
    else
        echo "Directory $1 already exists."
    fi
}

# Create directories using the function
create_dir_if_not_exists "data/amtDB"
create_dir_if_not_exists "data/mitogenomes_reich"
create_dir_if_not_exists "output"
create_dir_if_not_exists "src"
create_dir_if_not_exists "notebooks"

# Create Python package and scripts in src/
touch src/__init__.py
touch src/data_processing.py
touch src/sequence_analysis.py
touch src/utilities.py

# Initialize Jupyter notebook
touch notebooks/analysis_notebook.ipynb

# Create requirements.txt if it does not exist
if [ ! -f requirements.txt ]; then
    touch requirements.txt
    echo "Created requirements.txt file."
else
    echo "requirements.txt file already exists."
fi

# Create and activate virtual environment using Mamba
echo "Creating virtual environment with Mamba..."
mamba create -n $ENV_NAME 'python>=3.9' pandas biopython notebook -y
echo "To activate this environment, use:"
echo "     conda activate $ENV_NAME"

# Initialize Git repository if it does not exist
if [ ! -d ".git" ]; then
    echo "Initializing Git repository..."
    git init
    git add .
    git commit -m "Initial project setup"
    echo "Git repository initialized."
else
    echo "Git repository already exists."
fi

echo "Project setup complete. Your project is ready at: $(pwd)"
