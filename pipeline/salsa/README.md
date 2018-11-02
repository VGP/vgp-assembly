# Salsa2

Download [Salsa2.2](https://github.com/machinegun/SALSA/releases/tag/v2.2) in your $tools/salsa path

```
cd $tools
mkdir salsa2
wget https://github.com/machinegun/SALSA/archive/v2.2.tar.gz
tar -zxf v2.2.tar.gz
```

Build python environment: This needs to be done only once

```
mkdir -p $tools/conda/temp
cd $tools/conda/temp
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p $tools/conda -b
source $tools/conda/etc/profile.d/conda.sh
conda update conda
conda create -n salsa_env python=2.7.15 networkx==1.11
```

Test if salsa is working
```
source $tools/conda/etc/profile.d/conda.sh
conda activate salsa_env
python $tools/salsa2/SALSA-2.1/run_pipeline.py -h
```

