"""
CLI module for Molecular Docking Score Calculator.
"""
import argparse
import json
import sys

from .engine import (
    Atom, Molecule,
    calculate_vina_score, score_to_binding_affinity,
    binding_affinity_from_delta_g,
    lipinski_rule_of_5, extended_lipinski,
    estimate_admet, full_docking_analysis,
)


def _load_molecules(filepath: str):
    """Load ligand and receptor atoms from JSON."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    ligand_atoms = []
    for a in data.get('ligand', {}).get('atoms', []):
        ligand_atoms.append(Atom(
            element=a['element'], x=a['x'], y=a['y'], z=a['z'],
            atom_type=a.get('type', ''),
            is_donor=a.get('donor', False),
            is_acceptor=a.get('acceptor', False),
            is_hydrophobic=a.get('hydrophobic', False),
        ))
    
    receptor_atoms = []
    for a in data.get('receptor', {}).get('atoms', []):
        receptor_atoms.append(Atom(
            element=a['element'], x=a['x'], y=a['y'], z=a['z'],
            atom_type=a.get('type', ''),
            is_donor=a.get('donor', False),
            is_acceptor=a.get('acceptor', False),
            is_hydrophobic=a.get('hydrophobic', False),
        ))
    
    props = data.get('ligand', {}).get('properties', {})
    return ligand_atoms, receptor_atoms, props


def cmd_score(args):
    """Calculate docking score."""
    ligand_atoms, receptor_atoms, props = _load_molecules(args.input)
    num_rot = args.rotatable or props.get('rotatable_bonds', 0)
    
    score = calculate_vina_score(ligand_atoms, receptor_atoms, num_rot)
    affinity = score_to_binding_affinity(score.total)
    
    output = {
        'score_components': {
            'gauss1': score.gauss1,
            'gauss2': score.gauss2,
            'repulsion': score.repulsion,
            'hydrophobic': score.hydrophobic,
            'hbond': score.hbond,
            'torsional': score.torsional,
            'total': score.total,
        },
        'binding_affinity': {
            'delta_g': affinity.delta_g,
            'ki_nm': affinity.ki_nm,
            'classification': affinity.classification,
        },
    }
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print("Docking Score Analysis")
        print("=" * 50)
        print(f"Score Components:")
        print(f"  Gauss1 (attractive):  {score.gauss1:>10.4f}")
        print(f"  Gauss2 (repulsive):   {score.gauss2:>10.4f}")
        print(f"  Repulsion (steric):   {score.repulsion:>10.4f}")
        print(f"  Hydrophobic:          {score.hydrophobic:>10.4f}")
        print(f"  H-bond:               {score.hbond:>10.4f}")
        print(f"  Torsional:            {score.torsional:>10.4f}")
        print(f"  {'Total':>20}: {score.total:>10.4f} kcal/mol")
        print(f"\nBinding Affinity:")
        print(f"  ΔG:     {affinity.delta_g:.4f} kcal/mol")
        print(f"  Ki:     {affinity.ki_nm:.2f} nM")
        print(f"  Class:  {affinity.classification}")
    
    return 0


def cmd_lipinski(args):
    """Evaluate Lipinski Rule of 5."""
    result = lipinski_rule_of_5(args.mw, args.logp, args.hbd, args.hba)
    
    output = {
        'mw': result.mw, 'logp': result.logp,
        'hbd': result.hbd, 'hba': result.hba,
        'passes': result.passes,
        'violations': result.violations,
        'details': result.violation_details,
    }
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print("Lipinski Rule of 5")
        print("=" * 40)
        print(f"  MW:    {result.mw:.1f} {'✓' if result.mw <= 500 else '✗'} (≤500)")
        print(f"  LogP:  {result.logp:.2f} {'✓' if result.logp <= 5 else '✗'} (≤5)")
        print(f"  HBD:   {result.hbd} {'✓' if result.hbd <= 5 else '✗'} (≤5)")
        print(f"  HBA:   {result.hba} {'✓' if result.hba <= 10 else '✗'} (≤10)")
        print(f"\n  Result: {'PASS' if result.passes else 'FAIL'} ({result.violations} violations)")
        if result.violation_details:
            for v in result.violation_details:
                print(f"    - {v}")
    
    return 0


def cmd_admet(args):
    """Estimate ADMET properties."""
    result = estimate_admet(args.mw, args.logp, args.hbd, args.hba,
                             args.psa, args.rotatable)
    
    output = {
        'absorption': result.absorption,
        'distribution': result.distribution,
        'metabolism': result.metabolism,
        'excretion': result.excretion,
        'toxicity_risk': result.toxicity_risk,
        'oral_bioavailability': result.oral_bioavailability,
    }
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print("ADMET Property Estimation")
        print("=" * 40)
        print(f"  Absorption:           {result.absorption}")
        print(f"  Distribution:         {result.distribution}")
        print(f"  Metabolism:           {result.metabolism}")
        print(f"  Excretion:            {result.excretion}")
        print(f"  Toxicity Risk:        {result.toxicity_risk}")
        print(f"  Oral Bioavailability: {result.oral_bioavailability}")
    
    return 0


def cmd_affinity(args):
    """Calculate binding affinity from ΔG."""
    affinity = binding_affinity_from_delta_g(args.delta_g, args.temperature)
    
    output = {
        'delta_g': affinity.delta_g,
        'ki': affinity.ki,
        'ki_nm': affinity.ki_nm,
        'kd': affinity.kd,
        'classification': affinity.classification,
    }
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print("Binding Affinity")
        print("=" * 40)
        print(f"  ΔG:     {affinity.delta_g:.4f} kcal/mol")
        print(f"  Ki:     {affinity.ki:.4e} M")
        print(f"  Ki:     {affinity.ki_nm:.2f} nM")
        print(f"  Kd:     {affinity.kd:.4e} M")
        print(f"  Class:  {affinity.classification}")
    
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='molecular-docking-vina-scorer',
        description='Molecular Docking Score Calculator: Vina-like scoring, Lipinski, ADMET',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Score
    p_score = subparsers.add_parser('score', help='Calculate docking score')
    p_score.add_argument('-i', '--input', required=True, help='JSON with ligand/receptor atoms')
    p_score.add_argument('--rotatable', type=int, default=0, help='Number of rotatable bonds')
    p_score.add_argument('--json', action='store_true', help='JSON output')
    
    # Lipinski
    p_lip = subparsers.add_parser('lipinski', help='Lipinski Rule of 5')
    p_lip.add_argument('--mw', type=float, required=True, help='Molecular weight')
    p_lip.add_argument('--logp', type=float, required=True, help='LogP')
    p_lip.add_argument('--hbd', type=int, required=True, help='H-bond donors')
    p_lip.add_argument('--hba', type=int, required=True, help='H-bond acceptors')
    p_lip.add_argument('--json', action='store_true', help='JSON output')
    
    # ADMET
    p_admet = subparsers.add_parser('admet', help='ADMET estimation')
    p_admet.add_argument('--mw', type=float, required=True)
    p_admet.add_argument('--logp', type=float, required=True)
    p_admet.add_argument('--hbd', type=int, required=True)
    p_admet.add_argument('--hba', type=int, required=True)
    p_admet.add_argument('--psa', type=float, default=80.0)
    p_admet.add_argument('--rotatable', type=int, default=3)
    p_admet.add_argument('--json', action='store_true', help='JSON output')
    
    # Affinity
    p_aff = subparsers.add_parser('affinity', help='Binding affinity from ΔG')
    p_aff.add_argument('--delta-g', type=float, required=True, help='ΔG in kcal/mol')
    p_aff.add_argument('--temperature', type=float, default=298.15, help='Temperature (K)')
    p_aff.add_argument('--json', action='store_true', help='JSON output')
    
    args = parser.parse_args(argv)
    
    commands = {
        'score': cmd_score,
        'lipinski': cmd_lipinski,
        'admet': cmd_admet,
        'affinity': cmd_affinity,
    }
    
    return commands[args.command](args)


if __name__ == '__main__':
    sys.exit(main())
