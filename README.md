# Molecular Docking Score Calculator

AutoDock Vina-like scoring function, binding affinity prediction, Lipinski Rule of 5, and simplified ADMET estimation.

## Features

- **Vina-like Scoring**: Gauss1 (attractive), Gauss2 (repulsive), repulsion (steric clashes), hydrophobic contacts, H-bond geometry, torsional penalty
- **Binding Affinity**: ΔG → Ki conversion using Ki = e^(ΔG/RT)
- **Lipinski Rule of 5**: MW ≤500, LogP ≤5, HBD ≤5, HBA ≤10
- **Extended Rules**: Veber criteria (rotatable bonds, PSA)
- **ADMET Estimation**: Simplified absorption, distribution, metabolism, excretion, toxicity

## Quick Start

```bash
# Calculate docking score from atom coordinates
python cli.py score -i complex.json

# Lipinski Rule of 5
python cli.py lipinski --mw 350 --logp 2.5 --hbd 2 --hba 4

# ADMET estimation
python cli.py admet --mw 350 --logp 2.5 --hbd 2 --hba 4 --psa 80

# Binding affinity from ΔG
python cli.py affinity --delta-g -7.0
```

## Input Format (score command)

```json
{
  "ligand": {
    "atoms": [
      {"element": "C", "x": 0, "y": 0, "z": 0, "hydrophobic": true},
      {"element": "N", "x": 1.5, "y": 0, "z": 0, "donor": true}
    ],
    "properties": {"rotatable_bonds": 2}
  },
  "receptor": {
    "atoms": [
      {"element": "C", "x": 3, "y": 0, "z": 0, "hydrophobic": true},
      {"element": "O", "x": 3, "y": 2, "z": 0, "acceptor": true}
    ]
  }
}
```

## Python API

```python
from vina_docking_scorer import (
    Atom, calculate_vina_score, score_to_binding_affinity,
    lipinski_rule_of_5, estimate_admet,
)

# Score a docking pose
ligand = [Atom('C', 0, 0, 0, is_hydrophobic=True)]
receptor = [Atom('C', 3, 0, 0, is_hydrophobic=True)]
score = calculate_vina_score(ligand, receptor, num_rotatable_bonds=2)

# Binding affinity
affinity = score_to_binding_affinity(score.total)
print(f"Ki = {affinity.ki_nm:.1f} nM ({affinity.classification})")

# Drug-likeness
lip = lipinski_rule_of_5(mw=350, logp=2.5, hbd=2, hba=4)
print(f"Lipinski: {'PASS' if lip.passes else 'FAIL'}")
```

## Scoring Function

| Component | Weight | Description |
|-----------|--------|-------------|
| Gauss1 | -0.0356 | Attractive van der Waals (Gaussian at 3.0 Å) |
| Gauss2 | -0.00516 | Repulsive van der Waals (broader Gaussian) |
| Repulsion | 0.840 | Steric clashes (quadratic penalty < 2.0 Å) |
| Hydrophobic | -0.0351 | Hydrophobic contacts (1.5-4.0 Å ramp) |
| H-bond | -0.587 | H-bond geometry (optimal 2.8 Å, 180°) |
| Torsional | 0.0585 | Per rotatable bond penalty |

## License

MIT License.
