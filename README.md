# Molecular Docking Vina Scorer

> **Domain:** Computational Biology & AI Drug Discovery
> **Reference Guidelines & Standards:** `wwPDB, IUPAC & CLSI Computational Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**Molecular Docking Vina Scorer** is a computational platform implementing AutoDock Vina scoring potential (Gauss1/2, repulsion, H-bonding) for receptor binding prediction. It provides:

- **Vina-like Scoring Engine**: Deterministic calculation of docking scores using Gaussian-based potentials
- **Binding Affinity Prediction**: Convert docking scores to Ki/Kd values with classification
- **Lipinski Rule of 5**: Drug-likeness evaluation with violation details
- **ADMET Estimation**: Simplified absorption, distribution, metabolism, excretion, and toxicity prediction
- **Enterprise Audit System**: HMAC-SHA256 tamper-evident audit trail with PHI outbound guard

---

## ⚙️ Key Capabilities & Algorithmic Modules

- **Deterministic Calculation Engine**: Vina scoring function implementation with Gauss1/2, repulsion, hydrophobic, and hydrogen bonding terms
- **Binding Affinity Conversion**: ΔG to Ki/Kd calculation with classification (very_strong to very_weak)
- **Drug-Likeness Evaluation**: Lipinski Rule of 5 and extended Veber rules
- **ADMET Property Estimation**: Empirical rules for absorption, distribution, metabolism, excretion, and toxicity
- **Zero-PHI Outbound Guard**: Active regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers
- **Tamper-Evident Audit Trail**: Chained HMAC-SHA256 signed logs for every evaluation

---

## 🚀 Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/molecular-docking-vina-scorer.git
cd molecular-docking-vina-scorer

# Install dependencies
pip install -e .

# For development with tests
pip install -e ".[test]"
```

---

## 💻 CLI Quickstart & Usage

### 1. Calculate Docking Score
```bash
python cli.py score -i sample_payload.json
python cli.py score -i sample_payload.json --json
```

### 2. Evaluate Lipinski Rule of 5
```bash
python cli.py lipinski --mw 300 --logp 2.5 --hbd 2 --hba 4
python cli.py lipinski --mw 300 --logp 2.5 --hbd 2 --hba 4 --json
```

### 3. Estimate ADMET Properties
```bash
python cli.py admet --mw 300 --logp 2.5 --hbd 2 --hba 4 --psa 80 --rotatable 3
python cli.py admet --mw 300 --logp 2.5 --hbd 2 --hba 4 --json
```

### 4. Calculate Binding Affinity from ΔG
```bash
python cli.py affinity --delta-g -7.0
python cli.py affinity --delta-g -7.0 --temperature 298.15 --json
```

### 5. Enterprise Audit Commands
```bash
# Run task through enterprise supervisor
python cli.py audit --task-id TASK-001

# Query supervisory chat
python cli.py chat "Explain specifications"

# Verify audit trail integrity
python cli.py verify-audit
```

### Input Data Schema (for score command)

The score command expects a JSON file with the following structure:

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

---

## 🌐 REST API Server

Start the FastAPI server:

```bash
python -m uvicorn vina_docking_scorer.server:create_app --factory --host 0.0.0.0 --port 8000
```

Or using Docker:

```bash
docker-compose up
```

### API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/health` | Health check |
| POST | `/api/audit` | Process task through enterprise supervisor |
| POST | `/api/chat` | Query supervisory chat interface |

---

## 🛡️ Security Architecture

* **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, emails, DOBs, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Secure Key Management:** Audit signing key sourced from `AUDIT_SECRET_KEY` environment variable (no hardcoded secrets in production).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Run with coverage:

```bash
pytest -v --tb=short
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py 1000
```

---

## 🐳 Container Deployment

```bash
docker build -t molecular-docking-vina-scorer .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-secret-key molecular-docking-vina-scorer
```

Using docker-compose:

```bash
AUDIT_SECRET_KEY=your-secret-key docker-compose up
```

---

## 📁 Project Structure

```
molecular-docking-vina-scorer/
├── vina_docking_scorer/    # Core scoring engine
│   ├── engine.py           # Vina scoring, Lipinski, ADMET algorithms
│   ├── cli.py              # Command-line interface
│   ├── server.py           # FastAPI REST server
│   └── models.py           # Data model re-exports
├── agents/                 # Enterprise agent system
│   ├── base.py             # PHI guard, HMAC audit trail
│   ├── models.py           # Pydantic schemas
│   ├── supervisor.py       # Multi-agent orchestrator
│   ├── workers.py          # Specialized evaluation workers
│   ├── llm_factory.py      # LLM provider abstraction
│   └── ...
├── tests/                  # Test suite
├── web/                    # Operations console (HTML)
├── cli.py                  # Top-level CLI entry point
├── simulator.py            # Batch simulation
└── enrichment.py           # Domain enrichment features
```

---

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.
