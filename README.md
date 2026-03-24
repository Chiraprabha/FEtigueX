# FEtigueX

Phase field fatigue crack growth solver using Python API of Gmsh for geometry generation, FEniCSx for a staggered solver, and PyVista for visualization.

## Installation

Using conda/mamba:
```bash
conda env create -f environment.yml
conda activate fetiguex_env
pip install -e .
```

### IDE Setup

If your IDE shows import errors (e.g., missing `gmsh` or `yaml`), it is likely using your system's default Python instead of the Conda environment. 

To fix this locally:
1. **VS Code**: Press `Ctrl+Shift+P` (or `Cmd+Shift+P`), type **"Python: Select Interpreter"**, and select the `fetiguex_env` environment.
2. **PyCharm**: Go to `Settings > Project > Python Interpreter` and select the `fetiguex_env` environment.
