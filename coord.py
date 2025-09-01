import numpy as np
import re

# ====================== CONFIGURACIÓN ======================
gaussian_log = 'hoso.log'
wigner_file = 'wigner_gaussian.dat'
fit_coeff_file = 'fit_coefficients.dat'
nsamples =100 
bohr_to_angstrom = 0.529177
amu_to_au = 1822.888486

# ====================== LECTURA DE GEOMETRÍA ======================
coord_string = 'Ref.Geom.'
search_string = '-------------------'

coord = []
atoms = []

with open(gaussian_log) as f:
    coord_section = False
    count = 0
    for line in f:
        if coord_string in line:
            coord_section = True
            continue
        if coord_section:
            if count >= 4:
                if line.strip().startswith(search_string):
                    break
                if line.strip():
                    values = line.split()
                    atoms.append(int(values[1]))
                    coord.append([float(v) for v in values[2:]])
            count += 1

coord = np.array(coord)
natoms = len(atoms)
nmodes = natoms * 3 - 6


# ====================== LECTURA DE DESPLAZAMIENTOS ======================
displacements = []

for mode in range(1, nmodes + 1):
    vibration = []
    with open(gaussian_log) as f:
        inside_section = False
        countline = 0
        count = 0

        if (mode - 1) % 3 == 0:
            start, stop = 2, 5
        elif (mode - 2) % 3 == 0:
            start, stop = 5, 8
        else:
            start, stop = 8, 11

        linejump = ((mode - 1) // 3) * (natoms + 6)

        for line in f:
            if 'Atom  AN' in line:
                inside_section = True
                continue
            if inside_section:
                if countline >= linejump:
                    values = line.split()
                    if len(values) >= stop:
                        vibration.append([float(v) for v in values[start:stop]])
                        count += 1
                countline += 1
            if count >= natoms:
                break

    displacements.append(np.array(vibration))

# ====================== LECTURA DE a, sigma Y MASAS ======================
a_list = []
sigma_list = []
red_masses = []

with open(wigner_file) as f:
    for line in f:
        if line.strip().startswith("vib"):
            parts = line.split()
            a_list.append(float(parts[1]))
            sigma_list.append(float(parts[2]))

with open(fit_coeff_file) as f:
    for line in f:
        if line.strip().startswith("vib"):
            parts = line.split()
            red_masses.append(float(parts[4]))  # columna 5

a_list = np.array(a_list)
sigma_list = np.array(sigma_list)
red_masses = np.array(red_masses)

assert len(a_list) == nmodes and len(sigma_list) == nmodes
assert len(red_masses) == nmodes


# ====================== CONVERTIR a y sigma A·√amu -> Bohr/au -> Å ======================
# 1) dividir por sqrt(mu*amu_to_au) para pasar a unidades de displacement en a.u.
# 2) multiplicar por bohr_to_angstrom para Å
# ====================== AJUSTAR DESPLAZAMIENTOS DE GAUSSIAN ======================

for mode in range(nmodes):
    
    displacements[mode] /= np.sqrt(red_masses[mode])
    displacements[mode] *= bohr_to_angstrom  # ahora en Å   

# ====================== MUESTREO WIGNER ======================
def sample_wigner(coord, displacements, a_list, sigma_list, nsamples=5, seed=None):
    rng = np.random.default_rng(seed)
    samples = []

    for _ in range(nsamples):
        new_coords = np.array(coord, dtype=float)
        for mode in range(nmodes):
            q = rng.normal(a_list[mode], sigma_list[mode])
            new_coords += displacements[mode] * q
        samples.append(new_coords)
    return samples

samples = sample_wigner(coord, displacements, a_list, sigma_list, nsamples=nsamples)

# ====================== ESCRITURA ======================
atomic_dict = {}
with open('/home/macias/.scripts/Periodic_Table_mass.dat') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 4:
            atomic_dict[int(parts[0])] = parts[1]

symbols = [atomic_dict[Z] for Z in atoms]

with open('all_geometries.xyz', 'w') as fxyz:
    for i, geom in enumerate(samples, start=1):
        fxyz.write(f"{natoms}\n")
        fxyz.write(f"Wigner geometry {i}\n")
        for s, xyz in zip(symbols, geom):
            fxyz.write(f"{s:2s}  {xyz[0]:12.6f}  {xyz[1]:12.6f}  {xyz[2]:12.6f}\n")

print("generated geomtries in all_geometrias.xyz")
