"""Regression tests for public HTTP workflows and security boundaries."""
import json
import math

import pytest
from fastapi.testclient import TestClient

from agents.base import AuditTrail
from vina_docking_scorer.engine import Atom, calculate_vina_score, score_to_binding_affinity
from vina_docking_scorer.cli import main
from vina_docking_scorer.server import create_app


def _sample():
    return {
        "ligand": {
            "atoms": [{"element": "C", "x": 0.0, "y": 0.0, "z": 0.0, "hydrophobic": True}],
            "properties": {"rotatable_bonds": 5, "mw": 350, "logp": 2.5, "hbd": 2, "hba": 4},
        },
        "receptor": {
            "atoms": [{"element": "C", "x": 3.0, "y": 0.0, "z": 0.0, "hydrophobic": True}],
        },
    }


def test_http_application_and_browser_assets():
    app = create_app()
    assert app is not None
    client = TestClient(app)
    assert client.get("/health").status_code == 200
    assert client.get("/").status_code == 200
    js = client.get("/calculator.mjs")
    assert js.status_code == 200
    assert "export function analyze" in js.text


def test_http_scoring_and_input_validation():
    client = TestClient(create_app())
    response = client.post("/api/score", json=_sample())
    assert response.status_code == 200, response.text
    data = response.json()
    assert data["score"]["torsional"] == pytest.approx(0.2925)
    assert data["lipinski"]["passes"] is True
    assert data["admet"]["absorption"] == "high"

    bad = _sample()
    bad["ligand"]["atoms"][0]["x"] = "NaN"
    assert client.post("/api/score", json=bad).status_code == 422
    bad = _sample()
    bad["ligand"]["atoms"] = bad["ligand"]["atoms"] * 129
    assert client.post("/api/score", json=bad).status_code == 422
    bad = _sample()
    bad["temperature"] = 0
    assert client.post("/api/score", json=bad).status_code == 422


def test_audit_payload_blocked_with_client_error():
    client = TestClient(create_app())
    payload = {
        "task_id": "Patient MRN-123456", "target_identifier": "TEST",
        "primary_metric": 3.0,
    }
    assert client.post("/api/audit", json=payload).status_code == 400


def test_audit_signature_and_chain_tampering_detected():
    ledger = AuditTrail(secret_key="test-only-not-a-production-secret")
    ledger.log("actor", "tier", "event", {"status": "OK"})
    assert ledger.verify_integrity()
    ledger.get_trail()[0]["event_type"] = "external modification"
    assert ledger.verify_integrity()  # defensive snapshot

    ledger.logs[0]["event_type"] = "forged"
    assert not ledger.verify_integrity()
    ledger.logs[0]["event_type"] = "event"
    assert ledger.verify_integrity()
    ledger.logs[0]["current_hash"] = "0" * 64
    assert not ledger.verify_integrity()


def test_numeric_validation():
    with pytest.raises(ValueError):
        score_to_binding_affinity(-7.0, temperature=0)
    with pytest.raises(ValueError):
        score_to_binding_affinity(-7.0, temperature=math.nan)
    with pytest.raises(ValueError):
        score_to_binding_affinity(1e8)
    with pytest.raises(ValueError):
        calculate_vina_score([Atom("C", math.inf, 0, 0)], [Atom("C", 0, 0, 0)])
    with pytest.raises(ValueError):
        calculate_vina_score([], [], num_rotatable_bonds=-1)


def test_cli_zero_rotatable_override(tmp_path, capsys):
    source = tmp_path / "molecules.json"
    source.write_text(json.dumps(_sample()), encoding="utf-8")
    assert main(["score", "-i", str(source), "--rotatable", "0", "--json"]) == 0
    output = json.loads(capsys.readouterr().out)
    assert output["score_components"]["torsional"] == 0.0


def test_old_coordinator_is_importable_and_uses_supported_models():
    from agents.models import SystemTaskPayload
    from vina_docking_scorer.agents import MolecularDockingCoordinator
    coordinator = MolecularDockingCoordinator()
    report = coordinator.process(SystemTaskPayload(
        task_id="COMPAT-1", target_identifier="CHEM-1", primary_metric=2.0
    ))
    assert report["task_id"] == "COMPAT-1"
    assert report["overall_status"] == "ROUTINE"
    assert report["total_alerts"] == 0


def test_cli_rejects_non_docking_audit_payload(tmp_path):
    path = tmp_path / "audit-task.json"
    path.write_text(json.dumps({"task_id": "TASK-1", "primary_metric": 12}))
    with pytest.raises(ValueError, match="ligand and receptor"):
        main(["score", "-i", str(path)])


def test_audit_rejects_non_finite_measurements():
    client = TestClient(create_app())
    body = {"task_id": "A", "target_identifier": "B", "primary_metric": "NaN"}
    assert client.post("/api/audit", json=body).status_code == 422
