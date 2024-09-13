#!/bin/bash

function prepare_tsv() {
    local file_path="$1"

    # Check if the file path is provided
    if [ -z "$file_path" ]; then
        echo "Usage: $0 path_to_tsv_file"
        exit 1
    fi

    # Check if the file exists
    if [ ! -f "$file_path" ]; then
        echo "File not found!"
        exit 1
    fi

    # Create a temporary file to hold sorted output
    local temp_file="$(mktemp)"

    # Remove the first line (header), modify the first column, sort by the first column, and write to temp file
    {
        # Remove the first line (header)
        tail -n +2 "$file_path" |
        # Use awk to modify the first column (remove directory and extension) and keep tab-separated fields
        awk -F'\t' '{
            split($1, arr, "/");  # split by "/"
            filename = arr[length(arr)];  # last element is the filename with extension
            sub(/\.[^.]*$/, "", filename);  # remove the file extension
            $1 = filename;  # Keep only the filename without extension
            OFS="\t"; print $0;  # Preserve tabs between fields
        }' |
        # Replace all spaces with tabs (in case of some leftover spaces left)
        sed 's/ /	/g' |
        # Sort by the first column
        sort -t$'\t' -k1,1
    } > "$temp_file"

    # Add header back and overwrite the original file
    { head -n 1 "$file_path"; cat "$temp_file"; } > "$file_path"

    # Clean up
    rm "$temp_file"

    echo "File has been preprocessed: $file_path"
}

# Begin the script execution

# Call the function and pass the file path as an argument
prepare_tsv "$1"
