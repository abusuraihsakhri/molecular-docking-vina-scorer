"""
Molecular Docking Score Calculator: Vina-like scoring, binding affinity,
Lipinski Rule of 5, and ADMET estimation.
"""
__version__ = "3.0.0"

from .engine import (
    Atom, Molecule, ScoreComponents, BindingAffinity,
    LipinskiResult, ADMETEstimate, DockingResult,
    calculate_vina_score, score_to_binding_affinity,
    binding_affinity_from_delta_g,
    lipinski_rule_of_5, extended_lipinski,
    estimate_admet, full_docking_analysis,
    gauss1_score, gauss2_score, repulsion_score,
    hydrophobic_score, hbond_score,
)
