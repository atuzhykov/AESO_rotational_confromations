# MLIP Simulation Repository

This repository provides tools and scripts for simulating molecules using ASE calculators and GPAW.


### Інструкції з інсталяції MACE

If you do not have CUDA pre-installed, it is **recommended** to follow the conda installation process:

```sh
# Create a virtual environment and activate it
conda create --name mace_env
conda activate mace_env

# Install PyTorch
conda install pytorch torchvision torchaudio pytorch-cuda=11.6 -c pytorch -c nvidia

# (optional) Install MACE's dependencies from Conda as well
conda install numpy scipy matplotlib ase opt_einsum prettytable pandas e3nn

# Clone and install MACE (and all required packages)
git clone https://github.com/ACEsuit/mace.git
pip install ./mace
```



### Інструкції з інсталяції GPAW (Linux only)
- Крок 1: Створення Conda середовища

```sh
conda create -n gpaw python=3.8
conda activate gpaw
```
- Крок 2: Встановлення залежностей
```sh
conda install -c conda-forge libxc
conda install -c conda-forge openblas libblas
sudo apt-get update
sudo apt-get install gcc
```
- Крок 3: Налаштування змінних середовища
```sh
export C_INCLUDE_PATH=$(conda info --base)/envs/gpaw/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$(conda info --base)/envs/gpaw/include:$CPLUS_INCLUDE_PATH
export PYTHONPATH=/mnt/c/work/PINN-PoC/:$PYTHONPATH
```
- Крок 4: Встановлення GPAW
```sh
python3 -m pip install gpaw
```