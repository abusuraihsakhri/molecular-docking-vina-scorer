"""Compatibility adapters for the separate operational audit subsystem.

These classes use the supported agents.models and agents.supervisor APIs.
They do not perform molecular docking calculations.
"""
from typing import Any, Dict, List

from agents.models import SystemTaskPayload, AgentAlert
from agents.workers import (
    InvariantQCWorker,
    SafetyEscalationWorker,
    ProtocolConformanceWorker,
)
from agents.supervisor import SystemSupervisor


class BindingPocketGridAgent:
    def audit(self, payload: SystemTaskPayload) -> List[AgentAlert]:
        return InvariantQCWorker.evaluate(payload)


class VinaEnergyScoreAgent:
    def audit(self, payload: SystemTaskPayload) -> List[AgentAlert]:
        return SafetyEscalationWorker.evaluate(payload)


class TorsionalStrainPenaltyAgent:
    def audit(self, payload: SystemTaskPayload) -> List[AgentAlert]:
        return ProtocolConformanceWorker.evaluate(payload)


class MolecularDockingCoordinator:
    """Backward-compatible audit adapter, unrelated to pose scoring."""

    def __init__(self):
        self.supervisor = SystemSupervisor(model_provider="mock")
        self.execution_ledger: Dict[str, Dict[str, Any]] = {}

    def process(self, payload: SystemTaskPayload) -> Dict[str, Any]:
        if isinstance(payload, dict):
            payload = SystemTaskPayload.model_validate(payload)
        dossier = self.supervisor.process_task(payload)
        warnings = sum(a.urgency.value == "ELEVATED_RISK" for a in dossier.alerts)
        result = {
            "system": dossier.system_slug,
            "domain": dossier.domain,
            "task_id": dossier.task_id,
            "target_identifier": dossier.target_identifier,
            "overall_status": dossier.overall_urgency.value,
            "total_alerts": dossier.total_alerts,
            "critical_count": dossier.critical_alerts_count,
            "warning_count": warnings,
            "alerts": [a.to_dict() for a in dossier.alerts],
            "standard_specification": dossier.standard_reference,
            "consensus_summary": dossier.consensus_summary,
        }
        self.execution_ledger[payload.task_id] = result
        return result

    def query_supervisory_chat(self, query: str) -> str:
        return self.supervisor.query_supervisory_chat(query)
