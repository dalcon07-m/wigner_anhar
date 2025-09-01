import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")  # use non-GUI backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ------------------- INPUT -------------------
main_dir = os.getcwd()   # main working directory
output_file = "energies.dat"  # file to store energies per vibration
coef_file = "fit_coefficients.dat"  # file to store fitted coefficients
pdf_file = "vib_plots.pdf"  # PDF for plots
gaussian_log = "hoso.log"   # Gaussian log with reduced masses
# ---------------------------------------------

# Regex pattern to extract SCF energy from Gaussian output
scf_pattern = re.compile(r"SCF Done:\s+E\(.+?\)\s*=\s*([-\d\.]+)")

# --- Read reduced masses from Gaussian log ---
red_masses = []
if os.path.exists(gaussian_log):
    with open(gaussian_log, "r") as flog:
        for line in flog:
            if "Red. masses --" in line:
                tail = line.split("--", 1)[-1]  # extract numbers after '--'
                nums = re.findall(r"[-+]?\d*\.\d+|\d+", tail)
                red_masses.extend([float(n) for n in nums])
else:
    print(f"Warning: {gaussian_log} does not exist; no reduced masses added.")

# List vibration folders in numerical order
vibs = [d for d in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, d)) and d.startswith("vib")]
vibs.sort(key=lambda x: int(re.findall(r'\d+', x)[0]))

coef_lines = []  # store output coefficient lines

# --- Prepare PDF for plots ---
pdf_pages = PdfPages(pdf_file)
fig = plt.figure(figsize=(8.27, 11.69))  # A4 size
fig.suptitle("Cuadratic(k) + Cubic(kc) + Quartic (kq) Potential Fitting", fontsize=16, y=0.98)
plots_per_page = 4
plot_count = 0

rm_idx = 0  # index for reduced masses

for vib_folder in vibs:
    vib_path = os.path.join(main_dir, vib_folder)
    results = []

    # Read energies for each scaling factor in vibration folder
    scalings = [d for d in os.listdir(vib_path) if os.path.isdir(os.path.join(vib_path, d))]
    for scale in scalings:
        scale_dir = os.path.join(vib_path, scale)
        log_file = os.path.join(scale_dir, f"{vib_folder}_{scale}.log")
        if not os.path.exists(log_file):
            continue
        energy = None
        with open(log_file, "r") as f:
            for line in f:
                match = scf_pattern.search(line)  # extract SCF energy
                if match:
                    energy = float(match.group(1))
                    break
        if energy is not None:
            results.append((scale, energy))

    # Save energies.dat sorted by displacement
    try:
        results.sort(key=lambda x: float(x[0]))
    except ValueError:
        results.sort(key=lambda x: x[0])

    out_path = os.path.join(vib_path, output_file)
    with open(out_path, "w") as out:
        for scale, energy in results:
            out.write(f"{scale} {energy}\n")

    # Skip vibration if no data
    if not results:
        continue

    # Prepare data for polynomial fit
    data = np.array(results, dtype=float)
    x = data[:, 0]  # displacements
    y = data[:, 1]  # energies

    # 4th-degree polynomial fit: y = a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0
    coeffs = np.polyfit(x, y, 4)
    a4, a3, a2, a1, a0 = coeffs
    k  = 2 * a2    # quadratic coefficient
    kc = 6 * a3    # cubic coefficient
    kq = 24 * a4   # quartic coefficient

    # Assign reduced mass if available
    red_mass = red_masses[rm_idx] if rm_idx < len(red_masses) else None
    rm_idx += 1

    # Store line for coefficient file
    if red_mass is not None:
        coef_lines.append(f"{vib_folder} {k:.6f} {kc:.6f} {kq:.6f} {red_mass:.6f}")
    else:
        coef_lines.append(f"{vib_folder} {k:.6f} {kc:.6f} {kq:.6f} NA")

    # ---- Plot data + polynomial fit ----
    plot_count += 1
    ax = plt.subplot(2, 2, plot_count)
    ax.plot(x, y, 'o', label='Data')

    x_fit = np.linspace(min(x), max(x), 200)
    y_fit = a0 + a1*x_fit + a2*x_fit**2 + a3*x_fit**3 + a4*x_fit**4
    ax.plot(x_fit, y_fit, '-', label='Fit')

    
    # Calculate R²
    y_model = a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4  # values at original x
    ss_res = np.sum((y - y_model)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - (ss_res / ss_tot)
    
    ax.set_xlabel("Displacement (Angstrom)")
    ax.set_ylabel("Energy (A.U)")
    ax.set_title(vib_folder)

    # Show fit coefficients and reduced mass in the plot
    txt = f"k={k:.3f}\nkc={kc:.3f}\nkq={kq:.3f}\nR²={r2:.4f}"
    if red_mass is not None:
        txt += f"\nμ={red_mass:.3f} amu"
    ax.text(0.02, 0.95, txt,
            transform=ax.transAxes, fontsize=9, va='top', ha='left',
            bbox=dict(facecolor='white', alpha=0.7))

    ax.legend(loc='upper right', fontsize=8)

    # Save page when 4 subplots are filled
    if plot_count == plots_per_page:
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        pdf_pages.savefig(fig)
        plt.close(fig)
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.suptitle("Cuadratic(k) + Cubic(kc) + Quartic (kq) Potential Fitting", fontsize=16, y=0.98)
        plot_count = 0

# Save last incomplete page
if plot_count > 0:
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    pdf_pages.savefig(fig)
    plt.close(fig)

pdf_pages.close()

# Write coefficients with reduced masses to file
with open(coef_file, "w") as f:
    f.write("Vib k kc kq RedMass\n")
    for line in coef_lines:
        f.write(line + "\n")

print(f"PDF generated: {pdf_file}")
print(f"Coefficients saved in: {coef_file}")

