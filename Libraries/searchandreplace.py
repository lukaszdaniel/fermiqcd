import os
import shutil

# Get input from user
sin = input("Pattern to replace: ")
sout = input("Replace with: ")

# Iterate over all files in the current directory
for filename in os.listdir("."):
    if filename == "mdp_compatibility_macros.h":
        continue
    if not filename.endswith(".h"):
        continue

    print("Processing file:", filename)

    # Create a backup of the file
    backup_filename = filename + ".bak"
    shutil.copy(filename, backup_filename)

    try:
        with open(filename, "r") as file:
            lines = file.readlines()
    except Exception as e:
        print(f"Error opening {filename}: {e}")
        continue

    updated_lines = []
    modified = False

    for line in lines:
        if sin in line:
            print("<", line.strip())  # Print original line without extra newlines
            line = line.replace(sin, sout)
            print(">", line.strip())  # Print modified line
            modified = True
        updated_lines.append(line)

    # Only rewrite the file if changes were made
    if modified:
        with open(filename, "w") as file:
            file.writelines(updated_lines)

    print(f"Finished processing {filename}\n")

