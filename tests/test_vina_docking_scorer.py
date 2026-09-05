"""
Real tests for Molecular Docking Score Calculator.
Tests Vina-like scoring, binding affinity, Lipinski, and ADMET.
"""
import sys
import math
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from vina_docking_scorer.engine import (
    Atom, Molecule, ScoreComponents, BindingAffinity,
    LipinskiResult, ADMETEstimate, DockingResult,
    calculate_vina_score, score_to_binding_affinity,
    binding_affinity_from_delta_g,
    lipinski_rule_of_5, extended_lipinski,
    estimate_admet, full_docking_analysis,
    gauss1_score, gauss2_score, repulsion_score,
    hydrophobic_score, hbond_score,
    _distance,
    W_GAUSS1, W_GAUSS2, W_REPULSION, W_HYDROPHOBIC, W_HBOND, W_NUM_ROT,
)


# --- Helper functions ---

def _make_atom(element, x, y, z, donor=False, acceptor=False, hydrophobic=False):
    return Atom(element=element, x=x, y=y, z=z,
                is_donor=donor, is_acceptor=acceptor, is_hydrophobic=hydrophobic)


def _make_simple_ligand():
    """Create a simple ligand with a few atoms."""
    return [
        _make_atom('C', 0, 0, 0, hydrophobic=True),
        _make_atom('N', 1.5, 0, 0, donor=True),
        _make_atom('O', 0, 1.5, 0, acceptor=True),
        _make_atom('C', -1.5, 0, 0, hydrophobic=True),
    ]


def _make_simple_receptor():
    """Create a simple receptor binding site."""
    return [
        _make_atom('C', 3, 0, 0, hydrophobic=True),
        _make_atom('O', 3, 2, 0, acceptor=True),
        _make_atom('N', 3, -1, 1, donor=True),
        _make_atom('C', 4, 1, -1, hydrophobic=True),
    ]


# --- Test Distance ---

def test_distance_basic():
    """Test basic distance calculation."""
    a = _make_atom('C', 0, 0, 0)
    b = _make_atom('C', 3, 4, 0)
    assert abs(_distance(a, b) - 5.0) < 0.001


def test_distance_same_point():
    """Distance between same point should be 0."""
    a = _make_atom('C', 1, 2, 3)
    assert _distance(a, a) == 0.0


# --- Test Scoring Components ---

def test_gauss1_at_optimal():
    """Gauss1 should be maximal at optimal distance."""
    v = gauss1_score(3.0)
    assert v > gauss1_score(1.0)
    assert v > gauss1_score(10.0)
    assert abs(v - 1.0) < 0.001


def test_gauss1_symmetric():
    """Gauss1 should be symmetric around optimal."""
    v1 = gauss1_score(2.0)
    v2 = gauss1_score(4.0)
    assert abs(v1 - v2) < 0.001


def test_gauss2_at_optimal():
    """Gauss2 should be maximal at optimal distance."""
    v = gauss2_score(3.0)
    assert abs(v - 1.0) < 0.001


def test_gauss2_broader():
    """Gauss2 should be broader than Gauss1."""
    v1 = gauss1_score(5.0)
    v2 = gauss2_score(5.0)
    assert v2 > v1  # Gauss2 decays slower


def test_repulsion_no_clash():
    """No repulsion when atoms are far apart."""
    assert repulsion_score(5.0) == 0.0


def test_repulsion_clash():
    """Repulsion should be positive when atoms are too close."""
    v = repulsion_score(1.0)
    assert v > 0


def test_repulsion_at_cutoff():
    """At cutoff distance, repulsion should be 0."""
    assert repulsion_score(2.0) == 0.0


def test_hydrophobic_close():
    """Close hydrophobic atoms should score 1.0."""
    assert hydrophobic_score(1.0) == 1.0


def test_hydrophobic_far():
    """Far hydrophobic atoms should score 0.0."""
    assert hydrophobic_score(5.0) == 0.0


def test_hydrophobic_ramp():
    """Hydrophobic score should ramp linearly between good and bad."""
    v = hydrophobic_score(2.75)  # midpoint between 1.5 and 4.0
    assert abs(v - 0.5) < 0.01


def test_hbond_optimal():
    """H-bond at optimal distance and angle should score high."""
    v = hbond_score(2.8, 180.0)
    assert v > 0.9


def test_hbond_bad_angle():
    """H-bond with bad angle should score low."""
    v_good = hbond_score(2.8, 180.0)
    v_bad = hbond_score(2.8, 90.0)
    assert v_good > v_bad


# --- Test Vina Score ---

def test_vina_score_basic():
    """Basic Vina score calculation should return valid components."""
    ligand = _make_simple_ligand()
    receptor = _make_simple_receptor()
    score = calculate_vina_score(ligand, receptor, num_rotatable_bonds=2)
    
    assert isinstance(score, ScoreComponents)
    assert score.total != 0
    assert score.torsional == W_NUM_ROT * 2


def test_vina_score_no_atoms():
    """Empty atom lists should give zero score."""
    score = calculate_vina_score([], [])
    assert score.total == 0.0


def test_vina_score_torsional_penalty():
    """More rotatable bonds should increase (worsen) the score."""
    ligand = _make_simple_ligand()
    receptor = _make_simple_receptor()
    
    score_0 = calculate_vina_score(ligand, receptor, 0)
    score_5 = calculate_vina_score(ligand, receptor, 5)
    
    # Torsional penalty should make score_5 worse (more positive)
    assert score_5.torsional > score_0.torsional


def test_vina_score_close_atoms():
    """Very close atoms should have repulsion contribution."""
    ligand = [_make_atom('C', 0, 0, 0)]
    receptor = [_make_atom('C', 0.5, 0, 0)]  # very close
    score = calculate_vina_score(ligand, receptor)
    assert score.repulsion > 0


def test_vina_score_hbond():
    """Donor-acceptor pair should contribute H-bond score."""
    ligand = [_make_atom('N', 0, 0, 0, donor=True)]
    receptor = [_make_atom('O', 2.8, 0, 0, acceptor=True)]
    score = calculate_vina_score(ligand, receptor)
    assert score.hbond < 0  # favorable = negative


def test_vina_score_hydrophobic():
    """Hydrophobic pair should contribute hydrophobic score."""
    ligand = [_make_atom('C', 0, 0, 0, hydrophobic=True)]
    receptor = [_make_atom('C', 2.0, 0, 0, hydrophobic=True)]
    score = calculate_vina_score(ligand, receptor)
    assert score.hydrophobic < 0  # favorable = negative


# --- Test Binding Affinity ---

def test_binding_affinity_strong():
    """Very negative score should give strong binding."""
    affinity = score_to_binding_affinity(-10.0)
    assert affinity.classification in ('very_strong', 'strong')
    assert affinity.ki_nm < 1000


def test_binding_affinity_weak():
    """Slightly negative score should give weak binding."""
    affinity = score_to_binding_affinity(-1.0)
    assert affinity.classification in ('weak', 'very_weak')


def test_binding_affinity_positive_score():
    """Positive score should give very weak binding."""
    affinity = score_to_binding_affinity(1.0)
    assert affinity.classification == 'very_weak'


def test_binding_affinity_delta_g():
    """delta_g should equal the input score."""
    affinity = score_to_binding_affinity(-5.0)
    assert abs(affinity.delta_g - (-5.0)) < 0.001


def test_binding_affinity_from_delta_g():
    """binding_affinity_from_delta_g should work same as score_to_binding_affinity."""
    a1 = binding_affinity_from_delta_g(-7.0)
    a2 = score_to_binding_affinity(-7.0)
    assert abs(a1.ki_nm - a2.ki_nm) < 0.001


# --- Test Lipinski Rule of 5 ---

def test_lipinski_pass():
    """Drug-like molecule should pass Lipinski."""
    result = lipinski_rule_of_5(300, 2.5, 2, 4)
    assert result.passes
    assert result.violations == 0


def test_lipinski_fail_mw():
    """High MW should fail."""
    result = lipinski_rule_of_5(600, 2.5, 2, 4)
    assert not result.passes
    assert result.violations == 1


def test_lipinski_fail_logp():
    """High LogP should fail."""
    result = lipinski_rule_of_5(300, 6.0, 2, 4)
    assert not result.passes


def test_lipinski_fail_hbd():
    """Too many HBD should fail."""
    result = lipinski_rule_of_5(300, 2.5, 6, 4)
    assert not result.passes


def test_lipinski_fail_hba():
    """Too many HBA should fail."""
    result = lipinski_rule_of_5(300, 2.5, 2, 11)
    assert not result.passes


def test_lipinski_multiple_violations():
    """Multiple violations should be counted."""
    result = lipinski_rule_of_5(600, 6.0, 6, 11)
    assert result.violations == 4
    assert len(result.violation_details) == 4


def test_lipinski_boundary():
    """Exactly at boundary should pass."""
    result = lipinski_rule_of_5(500, 5.0, 5, 10)
    assert result.passes
    assert result.violations == 0


# --- Test Extended Lipinski ---

def test_extended_lipinski_pass():
    """Drug-like molecule should pass extended rules."""
    result = extended_lipinski(300, 2.5, 2, 4, psa=80, num_rotatable=5)
    assert result['overall_drug_like']


def test_extended_lipinski_fail_veber():
    """Too many rotatable bonds should fail Veber rules."""
    result = extended_lipinski(300, 2.5, 2, 4, psa=80, num_rotatable=15)
    assert not result['veber_passes']


# --- Test ADMET ---

def test_admet_drug_like():
    """Drug-like molecule should have favorable ADMET."""
    result = estimate_admet(300, 2.0, 2, 4, psa=80, num_rotatable=3)
    assert result.absorption in ('high', 'moderate')
    assert result.toxicity_risk in ('low', 'moderate')
    assert result.oral_bioavailability in ('high', 'moderate')


def test_admet_poor_drug():
    """Poor drug candidate should have unfavorable ADMET."""
    result = estimate_admet(700, 7.0, 8, 15, psa=200, num_rotatable=20)
    assert result.absorption == 'low'
    assert result.toxicity_risk in ('high', 'moderate')
    assert result.oral_bioavailability in ('low', 'very_low')


def test_admet_all_fields():
    """ADMET result should have all required fields."""
    result = estimate_admet(400, 3.0, 3, 6)
    assert result.absorption is not None
    assert result.distribution is not None
    assert result.metabolism is not None
    assert result.excretion is not None
    assert result.toxicity_risk is not None
    assert result.oral_bioavailability is not None


# --- Test Full Analysis ---

def test_full_docking_analysis():
    """Full analysis should return complete results."""
    ligand = _make_simple_ligand()
    receptor = _make_simple_receptor()
    result = full_docking_analysis(
        ligand, receptor,
        num_rotatable_bonds=2,
        mw=350, logp=2.5, hbd=2, hba=4,
    )
    assert isinstance(result, DockingResult)
    assert result.score.total != 0
    assert result.binding_affinity is not None
    assert result.lipinski is not None
    assert result.admet is not None


def test_full_docking_no_props():
    """Full analysis without molecular properties should skip Lipinski/ADMET."""
    ligand = _make_simple_ligand()
    receptor = _make_simple_receptor()
    result = full_docking_analysis(ligand, receptor)
    assert result.lipinski is None
    assert result.admet is None


# --- Test CLI ---

def test_cli_lipinski():
    """Test CLI lipinski command."""
    from vina_docking_scorer.cli import main
    assert main(['lipinski', '--mw', '300', '--logp', '2.5', '--hbd', '2', '--hba', '4']) == 0


def test_cli_lipinski_json():
    """Test CLI lipinski with JSON output."""
    from vina_docking_scorer.cli import main
    assert main(['lipinski', '--mw', '300', '--logp', '2.5', '--hbd', '2', '--hba', '4', '--json']) == 0


def test_cli_admet():
    """Test CLI admet command."""
    from vina_docking_scorer.cli import main
    assert main(['admet', '--mw', '300', '--logp', '2.5', '--hbd', '2', '--hba', '4']) == 0


def test_cli_affinity():
    """Test CLI affinity command."""
    from vina_docking_scorer.cli import main
    assert main(['affinity', '--delta-g', '-7.0']) == 0


def test_cli_score(tmp_path):
    """Test CLI score command."""
    data = {
        'ligand': {
            'atoms': [
                {'element': 'C', 'x': 0, 'y': 0, 'z': 0, 'hydrophobic': True},
                {'element': 'N', 'x': 1.5, 'y': 0, 'z': 0, 'donor': True},
            ],
            'properties': {'rotatable_bonds': 2},
        },
        'receptor': {
            'atoms': [
                {'element': 'C', 'x': 3, 'y': 0, 'z': 0, 'hydrophobic': True},
                {'element': 'O', 'x': 3, 'y': 2, 'z': 0, 'acceptor': True},
            ],
        },
    }
    f = tmp_path / 'docking.json'
    f.write_text(json.dumps(data))

    from vina_docking_scorer.cli import main
    assert main(['score', '-i', str(f)]) == 0


def test_cli_score_file_not_found():
    """Test CLI score with missing file raises FileNotFoundError."""
    from vina_docking_scorer.cli import main
    try:
        main(['score', '-i', '/nonexistent/path/file.json'])
        assert False, "Should have raised FileNotFoundError"
    except FileNotFoundError:
        pass


def test_cli_score_json_output(tmp_path):
    """Test CLI score with JSON output."""
    import io
    from contextlib import redirect_stdout

    data = {
        'ligand': {
            'atoms': [{'element': 'C', 'x': 0, 'y': 0, 'z': 0}],
            'properties': {},
        },
        'receptor': {
            'atoms': [{'element': 'C', 'x': 3, 'y': 0, 'z': 0}],
        },
    }
    f = tmp_path / 'docking.json'
    f.write_text(json.dumps(data))

    from vina_docking_scorer.cli import main

    captured = io.StringIO()
    with redirect_stdout(captured):
        result = main(['score', '-i', str(f), '--json'])

    assert result == 0
    output = json.loads(captured.getvalue())
    assert 'score_components' in output
    assert 'binding_affinity' in output
