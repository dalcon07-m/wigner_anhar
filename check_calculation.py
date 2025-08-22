import os

# ------------------- INPUT -------------------
main_dir = os.getcwd()   # Carpeta principal donde estás
# ---------------------------------------------

# Buscar todas las carpetas de vibraciones en la carpeta principal
vibs = [d for d in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, d)) and d.startswith("vib")]

for vib_folder in vibs:
    vib_path = os.path.join(main_dir, vib_folder)

    # Buscar todas las carpetas de scaling factor dentro de cada vib
    scalings = [d for d in os.listdir(vib_path) if os.path.isdir(os.path.join(vib_path, d))]

    for scale in scalings:
        scale_dir = os.path.join(vib_path, scale)

        # Buscar el archivo .log correspondiente
        log_file = os.path.join(scale_dir, f"{vib_folder}_{scale}.log")
        if not os.path.exists(log_file):
            print(f"[MISSING] {log_file} no encontrado")
            continue

        # Revisar si terminó correctamente
        with open(log_file, "r", errors="ignore") as f:
            content = f.read()

        if "Normal termination of Gaussian" in content:
            print(f"[OK] {vib_folder}_{scale} terminó correctamente")
        else:
            print(f"[FAIL] {vib_folder}_{scale} no terminó correctamente")

