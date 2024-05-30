#!/bin/bash

# Define the size limit (30 MB)
SIZE_LIMIT=30000

# Ensure .gitignore exists
touch .gitignore

# Find files larger than SIZE_LIMIT and add them to .gitignore
find . -type f -size +"${SIZE_LIMIT}"k | while read -r file; do
    # Normalize the file path to ensure consistency
    file=$(realpath --relative-to=. "$file")
    
    # Check if the file is already in .gitignore
    if ! grep -qxF "$file" .gitignore; then
        echo "$file" >> .gitignore
        echo "Added $file to .gitignore"
    else
        echo "$file is already in .gitignore"
    fi
    
    # Remove the file from the index if it is currently being tracked
    if [ -n "$file" ] && git ls-files --error-unmatch "$file" > /dev/null 2>&1; then
        git rm --cached "$file"
        echo "Removed $file from the index"
    fi
done

# Check if there are any changes to commit
if git diff --cached --quiet; then
    echo "No changes to commit."
else
    git commit -m "Add large files to .gitignore and stop tracking them"
    echo "Changes committed."
fi

echo "Script execution completed."
