for W in $(seq 24 3 60); do
    for i in $(seq 1 100); do
        sbatch <<EOF
#!/bin/bash
#SBATCH --account=m1266
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --time=14:00:00
#SBATCH --qos=regular
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --gpus-per-node=0
#SBATCH --job-name=W${W}_run${i}
module load julia
srun julia test_run.jl ${W} 1
EOF
    done
done