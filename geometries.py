import numpy as np  # import NumPy for numerical operations
import os  # import os for file and directory operations

# ====================== INPUT ======================
filename = 'hoso.log'  # Gaussian output file to read
scaling_factors = [-0.2, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20]  # scale factors for vibrational displacements

#____________________________________________________
#-----------------------------------------------------

coord_string = 'Ref.Geom.'  # marker for start of reference geometry in file
search_string = '-------------------'  # marker for end of geometry section

# ====================== FUNCTIONS ======================
def scale_and_save(matrix, mode_index, scale_factor, base_coord, atoms_symbols, outdir):
    """Scale displacement vectors and save as XYZ file in a given directory"""
    os.makedirs(outdir, exist_ok=True)  # create directory if it doesn't exist
    scaled_matrix = np.array(base_coord) + np.array(matrix) * scale_factor  # apply scaling to base coordinates
    outfile = os.path.join(outdir, f'vib{mode_index}_{scale_factor:.2f}.xyz')  # filename for this mode & scale
    with open(outfile, 'w') as f:  # open file for writing
        natoms = len(atoms_symbols)  # number of atoms
        f.write(f"{natoms}\n")  # write number of atoms
        f.write("Angstrom\n")  # write units
        for i in range(natoms):
            line = f"{atoms_symbols[i]} " + " ".join(f"{x:.6f}" for x in scaled_matrix[i])  # atom symbol + coordinates
            f.write(line + "\n")  # write line to file

# ====================== READ BASE GEOMETRY ======================
coord = []  # list to store coordinates
atoms = []  # list to store atomic numbers

with open(filename) as f:
    coord_section = False  # flag to detect geometry section
    count = 0
    for line in f:
        if coord_string in line:  # start of geometry section
            coord_section = True
            continue
        if coord_section:
            if count >= 4:  # skip Gaussian header lines
                if line.strip().startswith(search_string):  # end of geometry
                    break
                if line.strip():  # non-empty line
                    values = line.split()  # split line into elements
                    atoms.append(values[1])  # append atomic number
                    coord.append([float(v) for v in values[2:]])  # append coordinates
            count += 1

natoms = len(atoms)  # total number of atoms
nmodes = natoms * 3 - 6  # number of vibrational modes (non-linear molecule)

# ====================== ATOMIC SYMBOLS ======================
atomic_dict = {}  # dictionary to map atomic numbers to symbols
with open('/home/macias/.scripts/Periodic_Table.dat', 'r') as file:
    atomic_dict = {int(line.split()[0]): line.split()[1] for line in file}  # read periodic table

satoms = [atomic_dict[int(num)] for num in atoms]  # map atomic numbers to symbols

# ====================== READ VIBRATIONAL MODES AND SCALE ======================
for mode in range(1, nmodes + 1):
    vibration = []  # store displacement vectors for this mode
    with open(filename) as f:
        inside_section = False  # flag for mode section
        countline = 0  # track line number
        count = 0  # count atoms read for this mode

        # Determine which columns to read based on mode (Gaussian 3-column pattern)
        if (mode - 1) % 3 == 0:
            start, stop = 2, 5
        elif (mode - 2) % 3 == 0:
            start, stop = 5, 8
        else:
            start, stop = 8, 11

        # Calculate how many lines to skip to reach this mode
        linejump = ((mode - 1) // 3) * (natoms + 6)

        for line in f:
            if 'Atom  AN' in line:  # start of vibrational mode section
                inside_section = True
                continue
            if inside_section:
                if countline >= linejump:  # skip header lines
                    values = line.split()
                    if len(values) >= stop:  # ensure line has enough columns
                        vibration.append([float(v) for v in values[start:stop]])  # append displacements
                        count += 1
                countline += 1
            if count >= natoms:  # stop after reading all atoms
                break

    # Create directory for this mode
    mode_dir = f"vib{mode}"
    os.makedirs(mode_dir, exist_ok=True)

    # Scale and save each displacement in its subfolder
    for scale_factor in scaling_factors:
        scale_dir = os.path.join(mode_dir, f"{scale_factor:.2f}")
        scale_and_save(vibration, mode, scale_factor, coord, satoms, scale_dir)

