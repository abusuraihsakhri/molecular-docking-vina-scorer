"""FastAPI entry point for the calculator and the existing audit API.

The audit endpoints are implemented by agents.api. The scoring endpoint uses
the deterministic, simplified energy model in vina_docking_scorer.engine.
"""
from dataclasses import asdict
from pathlib import Path
from typing import List, Optional

from fastapi import HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

from agents.api import app
from .engine import Atom, full_docking_analysis


class AtomInput(BaseModel):
    element: str = Field(min_length=1, max_length=3)
    x: float = Field(allow_inf_nan=False)
    y: float = Field(allow_inf_nan=False)
    z: float = Field(allow_inf_nan=False)
    type: str = ""
    donor: bool = False
    acceptor: bool = False
    hydrophobic: bool = False


class LigandProperties(BaseModel):
    rotatable_bonds: int = Field(default=0, ge=0, le=500)
    mw: float = Field(default=0, ge=0, allow_inf_nan=False)
    logp: float = Field(default=0, allow_inf_nan=False)
    hbd: int = Field(default=0, ge=0)
    hba: int = Field(default=0, ge=0)
    psa: float = Field(default=80, ge=0, allow_inf_nan=False)


class MoleculeInput(BaseModel):
    atoms: List[AtomInput] = Field(min_length=1, max_length=128)
    properties: LigandProperties = Field(default_factory=LigandProperties)


class ScoreRequest(BaseModel):
    ligand: MoleculeInput
    receptor: MoleculeInput
    rotatable_bonds: Optional[int] = Field(default=None, ge=0, le=500)
    temperature: float = Field(default=298.15, gt=0, le=2000, allow_inf_nan=False)


def _atom(atom: AtomInput) -> Atom:
    return Atom(
        element=atom.element, x=atom.x, y=atom.y, z=atom.z,
        atom_type=atom.type, is_donor=atom.donor,
        is_acceptor=atom.acceptor, is_hydrophobic=atom.hydrophobic,
    )


@app.post("/api/score")
def api_score(request: ScoreRequest):
    props = request.ligand.properties
    rotatable = (
        request.rotatable_bonds if request.rotatable_bonds is not None
        else props.rotatable_bonds
    )
    try:
        result = full_docking_analysis(
            ligand_atoms=[_atom(a) for a in request.ligand.atoms],
            receptor_atoms=[_atom(a) for a in request.receptor.atoms],
            num_rotatable_bonds=rotatable, mw=props.mw, logp=props.logp,
            hbd=props.hbd, hba=props.hba, psa=props.psa,
            temperature=request.temperature,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    return asdict(result)


@app.get("/calculator.mjs", include_in_schema=False)
def browser_calculator():
    return FileResponse(
        Path(__file__).resolve().parent.parent / "web" / "calculator.mjs",
        media_type="text/javascript",
    )


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(Path(__file__).resolve().parent.parent / "web" / "index.html")


def create_app():
    """Factory-compatible application for existing Uvicorn/Docker commands."""
    return app
