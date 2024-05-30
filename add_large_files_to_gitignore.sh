#!/bin/bash

# Define the size limit (30 MB)
SIZE_LIMIT=30000

# Ensure .gitignore exists
touch .gitignore

# Find files larger than SIZE_LIMIT and add them to .gitignore
find . -type f -size +"${SIZE_LIMIT}"k | while read -r file; do
    # Check if the file is already in .gitignore
    if ! grep -qxF "$file" .gitignore; then
        echo >> .gitignore # Ensure new entry starts on a new line
        echo "$file" >> .gitignore
        echo "Added $file to .gitignore"
    else
        echo "$file is already in .gitignore"
    fi
done

echo "Script execution completed."
