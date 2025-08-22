import os

#####################################################################################
####################                 Input              #############################
#####################################################################################

base_dir = "."                 # Main folder containing the vibration folders
head_file = "hoso_head.com"    # Header file located in the main folder

#####################################################################################
####################   Function to create Gaussian .com files   #####################
#####################################################################################
def create_gaussian_input(head_file, xyz_file, output_file):
    """Create a Gaussian .com file by combining the header and the xyz coordinates"""
    # Read header
    with open(head_file, "r") as head:
        head_content = head.read().rstrip()  # remove trailing whitespace/newlines

    # Read xyz coordinates
    with open(xyz_file, "r") as xyz:
        lines = xyz.readlines()[2:]  # Skip first two lines of .xyz
        coords = "".join(line for line in lines if line.strip())  # remove empty lines

    # Write .com file
    with open(output_file, "w") as out:
        out.write(head_content + "\n")   # header + newline
        out.write(coords)                # coordinates (no extra newline between charge/multiplicity and first atom)
        out.write("\n")                  # only one final newline at the end


#####################################################################################
####################   Walk through all folders and create .com files   ############
#####################################################################################

for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith(".xyz"):
            xyz_path = os.path.join(root, file)
            com_path = os.path.splitext(xyz_path)[0] + ".com"  # Use the same name as .xyz
            create_gaussian_input(os.path.join(base_dir, head_file), xyz_path, com_path)
            print(f"Generated: {com_path}")

