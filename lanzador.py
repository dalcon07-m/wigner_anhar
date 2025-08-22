import os  # import os for file and directory operations

# ------------------- INPUT -------------------
main_dir = os.getcwd()   # main working directory
node = 1                 # number of nodes for PBS job
# --------------------------------------------

user = os.environ['USER']  # get current username from environment

# Find all vibration folders in the main directory
vibs = [d for d in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, d)) and d.startswith("vib")]

for vib_folder in vibs:
    vib_path = os.path.join(main_dir, vib_folder)  # full path to vibration folder

    # Find all scaling factor subfolders inside each vibration folder
    scalings = [d for d in os.listdir(vib_path) if os.path.isdir(os.path.join(vib_path, d))]

    for scale in scalings:
        scale_dir = os.path.join(vib_path, scale)  # full path to scaling folder

        # Look for the corresponding Gaussian .com input file
        com_file = os.path.join(scale_dir, f"{vib_folder}_{scale}.com")
        if not os.path.exists(com_file):
            print(f"Warning: {com_file} does not exist, skipping...")  # warn if input file is missing
            continue

        # Create a Torque/PBS job script in the same folder
        run_script = os.path.join(scale_dir, f"run_{vib_folder}_{scale}.sh")
        torque_content = f"""#!/bin/bash
#PBS -u {user}                     # set user
#PBS -N {vib_folder}_{scale}       # job name
#PBS -l nodes={node}:ppn=2         # request nodes and processors
#PBS -S /bin/bash                   # shell
#PBS -m be                          # email at begin and end
#PBS -r n                           # do not rerun on failure

FILE="{vib_folder}_{scale}"        # base file name
ScrDir="/scr/{user}/${{PBS_JOBID}}_$FILE"  # temporary scratch directory
Wdir="{scale_dir}"                 # working directory
. /soft/g09.c01/g09/bsd/g09.profile       # load Gaussian environment

mkdir -p $ScrDir                   # create scratch directory
cd $ScrDir                          # move to scratch
g09 < $Wdir/$FILE.com > $ScrDir/$FILE.log  # run Gaussian
cp *.log $Wdir                      # copy log files back
cp *.o $Wdir                         # copy .o files
cp *.e $Wdir                         # copy .e files
exit
"""
        # Write the Torque script to file
        with open(run_script, "w") as f:
            f.write(torque_content)

        # Submit job to the queue
        os.system(f"qsub {run_script}")
        print(f"Submitted job {vib_folder}_{scale} in {scale_dir}")  # print confirmation

