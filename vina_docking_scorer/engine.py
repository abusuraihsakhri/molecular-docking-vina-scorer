"""
Molecular Docking Score Calculator Engine
AutoDock Vina-like scoring, binding affinity prediction, Lipinski Rule of 5,
and simplified ADMET estimation.
"""
import math
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field


# --- Constants ---
R_GAS = 1.987e-3  # kcal/(mol·K)
T_DEFAULT = 298.15  # K (25°C)
RT = R_GAS * T_DEFAULT  # ~0.592 kcal/mol

# Vina scoring weights (from AutoDock Vina)
W_GAUSS1 = -0.0356
W_GAUSS2 = -0.00516
W_REPULSION = 0.840
W_HYDROPHOBIC = -0.0351
W_HBOND = -0.587
W_NUM_ROT = 0.0585  # torsional penalty per rotatable bond


# --- Data Models ---

@dataclass
class Atom:
    """Atom in 3D space."""
    element: str
    x: float
    y: float
    z: float
    atom_type: str = ''  # e.g., 'HD', 'OA', 'C', 'N', etc.
    is_donor: bool = False
    is_acceptor: bool = False
    is_hydrophobic: bool = False


@dataclass
class Molecule:
    """Simple molecule representation."""
    name: str
    atoms: List[Atom]
    num_rotatable_bonds: int = 0
    molecular_weight: float = 0.0
    logp: float = 0.0
    hbd: int = 0  # hydrogen bond donors
    hba: int = 0  # hydrogen bond acceptors
    num_heavy_atoms: int = 0


@dataclass
class ScoreComponents:
    """Individual scoring function components."""
    gauss1: float = 0.0      # Attractive van der Waals
    gauss2: float = 0.0      # Repulsive van der Waals
    repulsion: float = 0.0   # Steric clashes
    hydrophobic: float = 0.0 # Hydrophobic contacts
    hbond: float = 0.0       # Hydrogen bonding
    torsional: float = 0.0   # Torsional penalty
    total: float = 0.0


@dataclass
class BindingAffinity:
    """Binding affinity prediction."""
    delta_g: float           # kcal/mol
    ki: float                # Molar
    ki_nm: float             # Nanomolar
    kd: float                # Dissociation constant (M)
    classification: str      # strong, moderate, weak


@dataclass
class LipinskiResult:
    """Lipinski Rule of 5 evaluation."""
    mw: float
    logp: float
    hbd: int
    hba: int
    passes: bool
    violations: int
    violation_details: List[str]


@dataclass
class ADMETEstimate:
    """Simplified ADMET property estimation."""
    absorption: str          # high, moderate, low
    distribution: str
    metabolism: str
    excretion: str
    toxicity_risk: str       # low, moderate, high
    oral_bioavailability: str  # high, moderate, low, very_low
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class DockingResult:
    """Complete docking analysis result."""
    score: ScoreComponents
    binding_affinity: BindingAffinity
    lipinski: Optional[LipinskiResult] = None
    admet: Optional[ADMETEstimate] = None


# --- Scoring Functions ---

def _distance(a: Atom, b: Atom) -> float:
    """Euclidean distance between two atoms."""
    return math.sqrt(
        (a.x - b.x) ** 2 +
        (a.y - b.y) ** 2 +
        (a.z - b.z) ** 2
    )


def gauss1_score(d: float, optimal: float = 3.0) -> float:
    """Gauss1: attractive van der Waals term.
    Gaussian centered at optimal distance (~3.0 Å).
    """
    return math.exp(-((d - optimal) / 0.5) ** 2)


def gauss2_score(d: float, optimal: float = 3.0) -> float:
    """Gauss2: repulsive van der Waals term.
    Broader Gaussian for steric complementarity.
    """
    return math.exp(-((d - optimal) / 2.0) ** 2)


def repulsion_score(d: float, cutoff: float = 2.0) -> float:
    """Repulsion: steric clash penalty.
    Quadratic penalty when atoms are too close.
    """
    if d >= cutoff:
        return 0.0
    return (cutoff - d) ** 2


def hydrophobic_score(d: float, good: float = 1.5, bad: float = 4.0) -> float:
    """Hydrophobic contact score.
    Linear ramp: 1.0 at good distance, 0.0 at bad distance.
    """
    if d <= good:
        return 1.0
    elif d >= bad:
        return 0.0
    else:
        return 1.0 - (d - good) / (bad - good)


def hbond_score(d: float, angle: float = 180.0,
                distance_opt: float = 2.8, angle_opt: float = 180.0) -> float:
    """Hydrogen bond score based on distance and angle.
    Optimal: D-A distance ~2.8 Å, D-H-A angle ~180°.
    """
    dist_term = math.exp(-((d - distance_opt) / 0.7) ** 2)
    angle_term = math.cos(math.radians(angle - angle_opt)) ** 2
    return dist_term * angle_term


def calculate_vina_score(ligand_atoms: List[Atom],
                          receptor_atoms: List[Atom],
                          num_rotatable_bonds: int = 0) -> ScoreComponents:
    """Calculate AutoDock Vina-like score between ligand and receptor.
    
    Score = w1*Gauss1 + w2*Gauss2 + w3*Repulsion + w4*Hydrophobic + w5*HBond + w6*NumRot
    """
    if not isinstance(num_rotatable_bonds, int) or num_rotatable_bonds < 0:
        raise ValueError("num_rotatable_bonds must be a non-negative integer")
    if any(not all(math.isfinite(v) for v in (a.x, a.y, a.z))
           for a in (*ligand_atoms, *receptor_atoms)):
        raise ValueError("Atom coordinates must be finite")

    gauss1_total = 0.0
    gauss2_total = 0.0
    repulsion_total = 0.0
    hydrophobic_total = 0.0
    hbond_total = 0.0
    
    for lig_atom in ligand_atoms:
        for rec_atom in receptor_atoms:
            d = _distance(lig_atom, rec_atom)
            
            # Gauss terms (all atom pairs)
            gauss1_total += gauss1_score(d)
            gauss2_total += gauss2_score(d)
            
            # Repulsion (close contacts)
            repulsion_total += repulsion_score(d)
            
            # Hydrophobic (both hydrophobic)
            if lig_atom.is_hydrophobic and rec_atom.is_hydrophobic:
                hydrophobic_total += hydrophobic_score(d)
            
            # H-bond (donor-acceptor pairs)
            if (lig_atom.is_donor and rec_atom.is_acceptor) or \
               (lig_atom.is_acceptor and rec_atom.is_donor):
                hbond_total += hbond_score(d)
    
    torsional = num_rotatable_bonds
    
    total = (W_GAUSS1 * gauss1_total +
             W_GAUSS2 * gauss2_total +
             W_REPULSION * repulsion_total +
             W_HYDROPHOBIC * hydrophobic_total +
             W_HBOND * hbond_total +
             W_NUM_ROT * torsional)
    
    return ScoreComponents(
        gauss1=round(W_GAUSS1 * gauss1_total, 4),
        gauss2=round(W_GAUSS2 * gauss2_total, 4),
        repulsion=round(W_REPULSION * repulsion_total, 4),
        hydrophobic=round(W_HYDROPHOBIC * hydrophobic_total, 4),
        hbond=round(W_HBOND * hbond_total, 4),
        torsional=round(W_NUM_ROT * torsional, 4),
        total=round(total, 4),
    )


# --- Binding Affinity ---

def score_to_binding_affinity(score: float, temperature: float = T_DEFAULT) -> BindingAffinity:
    """Convert docking score to binding affinity.
    
    Score (kcal/mol) ≈ ΔG_binding
    Ki = e^(ΔG / RT)
    """
    if not math.isfinite(score) or not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("score must be finite and temperature must be positive and finite")
    rt = R_GAS * temperature
    delta_g = score  # kcal/mol, meaningful for calibrated free energies only
    exponent = delta_g / rt
    if not -700 <= exponent <= 680:
        raise ValueError("Energy/temperature combination is outside the numeric range")
    ki = math.exp(exponent)
    
    ki_nm = ki * 1e9  # convert to nanomolar
    kd = ki  # simplified: Kd ≈ Ki
    
    if ki_nm < 10:
        classification = 'very_strong'
    elif ki_nm < 100:
        classification = 'strong'
    elif ki_nm < 1000:
        classification = 'moderate'
    elif ki_nm < 10000:
        classification = 'weak'
    else:
        classification = 'very_weak'
    
    return BindingAffinity(
        delta_g=round(delta_g, 4),
        ki=ki,
        ki_nm=round(ki_nm, 4),
        kd=kd,
        classification=classification,
    )


def binding_affinity_from_delta_g(delta_g: float, temperature: float = T_DEFAULT) -> BindingAffinity:
    """Calculate binding affinity from free energy change."""
    return score_to_binding_affinity(delta_g, temperature)


# --- Lipinski Rule of 5 ---

def lipinski_rule_of_5(mw: float, logp: float, hbd: int, hba: int) -> LipinskiResult:
    """Evaluate Lipinski's Rule of 5 for drug-likeness.
    
    Criteria:
    - Molecular weight ≤ 500
    - LogP ≤ 5
    - H-bond donors ≤ 5
    - H-bond acceptors ≤ 10
    
    Returns pass/fail and violation details.
    """
    violations = []
    
    if mw > 500:
        violations.append(f"MW {mw:.1f} > 500")
    if logp > 5:
        violations.append(f"LogP {logp:.2f} > 5")
    if hbd > 5:
        violations.append(f"HBD {hbd} > 5")
    if hba > 10:
        violations.append(f"HBA {hba} > 10")
    
    return LipinskiResult(
        mw=mw,
        logp=logp,
        hbd=hbd,
        hba=hba,
        passes=len(violations) == 0,
        violations=len(violations),
        violation_details=violations,
    )


def extended_lipinski(mw: float, logp: float, hbd: int, hba: int,
                       psa: float = 0.0, num_rotatable: int = 0) -> Dict[str, Any]:
    """Extended drug-likeness rules beyond Lipinski.
    
    Veber rules:
    - Rotatable bonds ≤ 10
    - PSA ≤ 140 Å²
    """
    lipinski = lipinski_rule_of_5(mw, logp, hbd, hba)
    
    veber_violations = []
    if num_rotatable > 10:
        veber_violations.append(f"Rotatable bonds {num_rotatable} > 10")
    if psa > 140:
        veber_violations.append(f"PSA {psa:.1f} > 140")
    
    return {
        'lipinski': lipinski,
        'veber_passes': len(veber_violations) == 0,
        'veber_violations': veber_violations,
        'overall_drug_like': lipinski.passes and len(veber_violations) == 0,
    }


# --- Simplified ADMET Estimation ---

def estimate_admet(mw: float, logp: float, hbd: int, hba: int,
                    psa: float = 80.0, num_rotatable: int = 3) -> ADMETEstimate:
    """Simplified ADMET property estimation based on molecular descriptors.
    
    Uses empirical rules based on molecular properties.
    This is a simplified model for educational purposes.
    """
    # Absorption (based on LogP, PSA, MW)
    if logp < 0 or logp > 5 or psa > 140 or mw > 500:
        absorption = 'low'
    elif 1 <= logp <= 3 and psa < 100:
        absorption = 'high'
    else:
        absorption = 'moderate'
    
    # Distribution (based on LogP and protein binding likelihood)
    if logp < 0:
        distribution = 'low'
    elif logp > 3:
        distribution = 'high'
    else:
        distribution = 'moderate'
    
    # Metabolism (based on MW and functional groups)
    if mw > 500 or hba > 8:
        metabolism = 'high'  # more sites for metabolism
    elif mw < 300 and hba < 4:
        metabolism = 'low'
    else:
        metabolism = 'moderate'
    
    # Excretion (related to MW and polarity)
    if mw > 500 or psa > 120:
        excretion = 'low'  # harder to excrete large molecules
    elif mw < 300 and psa < 80:
        excretion = 'high'
    else:
        excretion = 'moderate'
    
    # Toxicity risk (simplified)
    toxicity_factors = 0
    if logp > 5:
        toxicity_factors += 1
    if mw > 500:
        toxicity_factors += 1
    if hba > 10:
        toxicity_factors += 1
    if psa < 20:
        toxicity_factors += 1  # very lipophilic
    
    if toxicity_factors >= 3:
        toxicity_risk = 'high'
    elif toxicity_factors >= 1:
        toxicity_risk = 'moderate'
    else:
        toxicity_risk = 'low'
    
    # Oral bioavailability (composite)
    bio_score = 0
    if 200 <= mw <= 500:
        bio_score += 1
    if 0 <= logp <= 3:
        bio_score += 1
    if hbd <= 5:
        bio_score += 1
    if hba <= 10:
        bio_score += 1
    if psa <= 140:
        bio_score += 1
    if num_rotatable <= 10:
        bio_score += 1
    
    if bio_score >= 5:
        oral_bioavailability = 'high'
    elif bio_score >= 3:
        oral_bioavailability = 'moderate'
    elif bio_score >= 1:
        oral_bioavailability = 'low'
    else:
        oral_bioavailability = 'very_low'
    
    return ADMETEstimate(
        absorption=absorption,
        distribution=distribution,
        metabolism=metabolism,
        excretion=excretion,
        toxicity_risk=toxicity_risk,
        oral_bioavailability=oral_bioavailability,
        details={
            'mw': mw, 'logp': logp, 'hbd': hbd, 'hba': hba,
            'psa': psa, 'num_rotatable': num_rotatable,
            'bio_score': bio_score,
        },
    )


# --- Full Docking Analysis ---

def full_docking_analysis(ligand_atoms: List[Atom],
                           receptor_atoms: List[Atom],
                           num_rotatable_bonds: int = 0,
                           mw: float = 0.0, logp: float = 0.0,
                           hbd: int = 0, hba: int = 0,
                           psa: float = 80.0,
                           temperature: float = T_DEFAULT) -> DockingResult:
    """Run complete docking analysis.
    
    1. Calculate Vina-like score
    2. Predict binding affinity
    3. Evaluate Lipinski Rule of 5
    4. Estimate ADMET properties
    """
    score = calculate_vina_score(ligand_atoms, receptor_atoms, num_rotatable_bonds)
    affinity = score_to_binding_affinity(score.total, temperature)
    
    lipinski = None
    admet = None
    if mw > 0:
        lipinski = lipinski_rule_of_5(mw, logp, hbd, hba)
        admet = estimate_admet(mw, logp, hbd, hba, psa, num_rotatable_bonds)
    
    return DockingResult(
        score=score,
        binding_affinity=affinity,
        lipinski=lipinski,
        admet=admet,
    )
