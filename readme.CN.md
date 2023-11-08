<a href="./readme.CN.md">简体中文</a> | <a href="./readme.md">English</a>

![License](https://img.shields.io/badge/license-MIT-yellowgreen)  ![Language](https://img.shields.io/badge/language-python-blue)

# ProtIon

蛋白和离子互作的分子动力学模拟

## 安装

```bash
# git clone本repo, 之后进入本repo目录后运行以下指令

conda install -yc conda-forge mamba
mamba install -yc conda-forge rdkit openmm openmmforcefields pdbfixer mdtraj openff-toolkit numpy nglview biopython
pip install .
```

## 使用
### 加入离子并运行分子动力学模拟

你可以查看 ``demo/run_test.py``

```python
from protion import ProtIon

input_pdb = f"demo/inputs/folding_result.pdb"
outdir = f"demo/outputs"
# 给出目标离子以及对应浓度(单位：mol / L)
# 目标离子可以只给元素名、SMILES；或是给出"{元素名/SMILES}::{最长3个字符}"的命名，即中间用"::"隔开
# 例如下面
ions_concentration = {
    "Ca": 8e-3, # 目标离子只给元素名称
    "Na": 1.6e-2,
    "Mg": 5e-2,
    "C(=O)(O)[O-]::BCR": 1.6e-2, # 目标离子的SMILES为C(=O)(O)[O-], 命名为BCR
    "Cl": 8e-3 * 2 + 5e-2 * 2,
}
pi = ProtIon(outdir)
modeller, forcefield = pi.preprocess(
    input_pdb = input_pdb, 
    ions_concentratsion = ions_concentration,
    ph = 7.
)
pi.run_md(
    modeller,
    forcefield,
    step_time: float = 4,
    total_time: float = 5e6,
    skip_time: int = 1e6,
    save_interval_time: int = 1e3,
)
```
### 离子富集能力打分
```python
# 假设你按照上述方式增加离子并运行完了分子动力学模拟
# 我们现在对Ca离子浓度进行残基富集能力打分
from protion import ion_enrichment_score
import mdtraj as md

traj = md.load(
    "demo/outputs/trajectory_nowat.xtc", 
    top = "demo/outputs/_result_nowat.pdb"
)

Ca_enrichment_score = ion_enrichment_score(
    traj: md.core.trajectory.Trajectory,
    ion_resname: str,
    distance_cutoff: float = 0.4, # nm
    scheme: str = "closest",
)

```