import os
import numpy as np
import matplotlib
matplotlib.use("Agg")  # use non-GUI backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit  # for robust Gaussian fitting

# ------------------- INPUTS -------------------
fit_coeff_file = "fit_coefficients.dat"  # file with k, kc, kq, reduced masses
pdf_file = "wave_wigner.pdf"            # output PDF for plots
gauss_file="wigner_gaussian.dat"        # output Gaussian fit parameters
wigner_dir = "wigner_gaussian"          # directory for storing outputs
amu2emass = 1822.89                     # atomic mass unit to electron mass
hbar = 1.0
L =10.0                                  # spatial grid half-length
Np = 500                                # number of points in grid
plots_per_page = 9                       # subplots per PDF page
# ---------------------------------------------

# ------------------- Parse fitting coefficients -------------------
modes = []
with open(fit_coeff_file) as f:
    lines = f.readlines()[1:]  # skip header
for line in lines:
    parts = line.strip().split()
    if len(parts) != 5:
        continue
    vib_id = parts[0]
    k, kc, kq, mr = map(float, parts[1:])
    modes.append({"vib": vib_id, "k": k, "kc": kc, "kq": kq, "mr": mr})

if not modes:
    print("No modes found in fit_coefficients.dat")
    exit()

# ------------------- Functions -------------------
def build_hamiltonian(k, kc, kq, mr_emass, xgrid):
    """Build 1D Hamiltonian matrix with anharmonic potential"""
    n = xgrid.size
    dx = xgrid[1]-xgrid[0]
    V = np.diag(0.5*k*xgrid**2 + (1/6)*kc*xgrid**3 + (1/24)*kq*xgrid**4)
    main = np.full(n, -2.0)
    off = np.ones(n-1)
    D2 = (np.diag(main)+np.diag(off,1)+np.diag(off,-1))/(dx*dx)
    T = - (hbar**2)/(2*mr_emass)*D2
    return T+V

def normalize_wave(psi, dx):
    """Normalize wavefunction"""
    return psi/np.sqrt(np.trapz(np.abs(psi)**2, dx=dx))

def compute_ground_state(k, kc, kq, mr, L=L, Np=Np):
    """Compute ground state wavefunction for given potential"""
    mr_emass = mr*(amu2emass)
    xgrid = np.linspace(-L, L, Np)
    H = build_hamiltonian(k, kc, kq, mr_emass, xgrid)
    evals, evecs = np.linalg.eigh(H)
    psi0 = normalize_wave(evecs[:, np.argmin(evals)], xgrid[1]-xgrid[0])
    return xgrid, psi0

def wigner_p0(psi, dx):
    """Compute Wigner function at P=0 (using convolution)"""
    psi_rev = psi[::-1]
    W = np.convolve(psi, psi_rev, mode='same')*dx/(np.pi*hbar)
    W_norm= W/np.sum(W*dx)
    return W_norm  # normalized

def gaussian(x, A, x0, sigma):
    """1D Gaussian function"""
    return A * np.exp(-(x-x0)**2/(2*sigma**2))

def fit_gaussian(x, y):
    """Fit Gaussian to data using non-linear least squares"""
    #y_norm = y/np.max(y)                 # normalize data
    p0 = [np.max(y),x[np.argmax(y)], np.std(x)]  # initial guess: [A,center, width]
    
    lower_bounds = [0, -np.inf, 1e-12]  # sigma mínimo pequeño
    upper_bounds = [np.inf, np.inf, np.inf]
    
    
    
    try:
        popt, _ = curve_fit(gaussian, x, y, p0=p0, bounds=(lower_bounds, upper_bounds))
        A, x0, sigma = popt
    except RuntimeError:
        # fallback if fit fails
        x0 = x[np.argmax(y)]
        sigma = np.sqrt(np.sum(y*(x-x0)**2)/np.sum(y))
        sigma = max(sigma, 1e-12)
        A=np.max(y)
    W_fit = gaussian(x, A, x0, sigma)
    
    # calculate R²
    ss_res = np.sum((y - W_fit)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot
    
    return W_fit, x0, sigma, r2, A

# ------------------- PDF setup -------------------
pdf_pages = PdfPages(pdf_file)
plot_count = 0
fig = plt.figure(figsize=(11,8))

# ------------------- Output file header -------------------
with open(gauss_file, "w") as gf:
    gf.write("# vib   a   sigma; (Gaussian fit: W(x) = A*exp(-(x-a)^2 / (2*sigma^2)) !! A doesn't matter for sample   ) \n")

# ------------------- MAIN LOOP -------------------
for m in modes:
    vib, k, kc, kq, mr = m["vib"], m["k"], m["kc"], m["kq"], m["mr"]

    # --- Compute ground state wavefunction (anharmonic) ---
    xgrid, psi0 = compute_ground_state(k, kc, kq, mr)
    psi_abs = np.abs(psi0)
    psi_abs_2=psi_abs**2
    # --- Compute Wigner function (anharmonic) ---
    x_wig, psi_full = compute_ground_state(k, kc, kq, mr)
    dx_wig = x_wig[1]-x_wig[0]
    W_anh = wigner_p0(psi_full, dx_wig)
     
    # --- Harmonic approximation ---
    k_harm = k
    x_harm, psi_h = compute_ground_state(k_harm, 0, 0, mr)
    psi_h_abs = np.abs(psi_h)
    W_harm = wigner_p0(psi_h, dx_wig)

    # --- Fit Gaussian robustly ---
    W_fit, x0, sigma, r2 , A = fit_gaussian(x_wig, W_anh)

    # --- Save Gaussian parameters ---
    with open(gauss_file, "a") as gf:
        gf.write(f"{vib} {x0:.12e}  {sigma:.12e}\n")

    # --- Plot |Psi0(x)| anharmonic vs harmonic ---
    plot_count +=1
    plt.subplot(3,3,plot_count)
    plt.plot(xgrid, psi_abs, 'o-', label='|Psi0(x)| Anharmonic')
    plt.plot(x_harm, psi_h_abs, '-', label='|Psi0(x)| Harmonic')
    plt.xlabel("x (Angstroms)")
    plt.ylabel("|Psi0(x)|")
    plt.title(f"{vib} Wavefunction")
    plt.legend(loc='upper right', fontsize=8)
    plt.xlim(-2,2)
    if plot_count == plots_per_page:
        plt.tight_layout()
        pdf_pages.savefig(fig)
        plt.close(fig)
        fig = plt.figure(figsize=(11,8))
        plot_count=0

    # --- Plot Wigner harmonic vs anharmonic ---
    plot_count +=1
    plt.subplot(3,3,plot_count)
    plt.plot(xgrid, W_harm, '-', label='Harmonic Wigner')
    plt.plot(xgrid, W_anh, '-', label='Anharmonic Wigner')
    plt.plot(xgrid, psi_abs_2, 'o-', label='|Psi0(x)|²', markersize=3 )
    plt.xlabel("x (Angstroms)")
    plt.ylabel("W(x, P=0)")
    plt.title(f"{vib} Wigner Harm vs Anharm")
    plt.legend(loc='upper right', fontsize=8)
    plt.xlim(-3,3)
    if plot_count == plots_per_page:
        plt.tight_layout()
        pdf_pages.savefig(fig)
        plt.close(fig)
        fig = plt.figure(figsize=(11,8))
        plot_count=0

    # --- Plot Wigner original vs Gaussian ---
    plot_count +=1
    plt.subplot(3,3,plot_count)
    plt.plot(xgrid, W_anh, 'o', label='Original Wigner')
    plt.plot(xgrid, W_fit, '-', label=f'Gaussian Fit, R²={r2:.4f}')
    plt.xlabel("x (Angstroms)")
    plt.ylabel("W(x, P=0)")
    plt.title(f"{vib} Wigner Fit")
    plt.legend(loc='upper right', fontsize=8)
    plt.xlim(-2,2)
    if plot_count == plots_per_page:
        plt.tight_layout()
        pdf_pages.savefig(fig)
        plt.close(fig)
        fig = plt.figure(figsize=(11,8))
        plot_count=0
# --- Save remaining plots ---
if plot_count > 0:
    plt.tight_layout()
    pdf_pages.savefig(fig)
    plt.close(fig)
    
pdf_pages.close() 
print(f"PDF generated: {pdf_file}")
print(f"Gaussian parameters saved in {gauss_file}/")
print(f"Ahora solo falta conectar con NX :)")
