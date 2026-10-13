# Molecular Docking Vina Scorer

A deterministic, **Vina-like pose scoring demonstrator** with a Python command-line interface, a FastAPI service, and an entirely client-side browser calculator.

The browser source is in [web/](web/). GitHub Pages deployment is defined in [pages.yml](.github/workflows/pages.yml); a public URL should only be added once deployment is enabled and verified.

## Scientific scope and limitations

Given **already positioned** ligand and receptor atoms (Cartesian Å coordinates), the project sums simplified Gaussian, steric repulsion, hydrophobic, and donor–acceptor interaction terms, plus a rotatable-bond penalty. It also provides Lipinski/Veber screening and descriptor-based ADMET heuristics.

**It does not perform AutoDock Vina docking**, conformer search, pose optimization, receptor preparation, atom typing based on chemistry, or validated affinity prediction. Its potential uses fixed distance centers rather than full Vina atom-type-specific surface distances. The apparent Ki/Kd values are illustrative applications of `exp(ΔG/RT)` to the heuristic score—not experimentally calibrated binding constants. Descriptor-based ADMET labels are not validated clinical or toxicological predictions. Use only for education, testing, and exploratory analysis.

## Browser calculator

Run locally with an HTTP server (ES modules should not be opened using `file://`):

```bash
python -m http.server 8080 --directory web
```

Open `http://localhost:8080/`. Paste a JSON structure containing ligand/receptor atoms, optionally open a local JSON file, calculate the score, inspect the components and estimated descriptors, and download the result as JSON. Each molecule is limited to **128 atoms**. All calculator operations run locally in JavaScript; no uploads or API keys are needed. The separate supervisory audit form requires the Python server and is unavailable on static GitHub Pages.

**Browser requirements:** modern browsers supporting JavaScript ES modules. Browser calculations are implemented in JavaScript rather than Pyodide, avoiding large Python/WASM downloads and native chemistry dependencies.

## Python install and CLI

Requires Python 3.9+:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
python cli.py score -i sample_docking.json --json
python cli.py score -i sample_docking.json --rotatable 0
python cli.py lipinski --mw 300 --logp 2.5 --hbd 2 --hba 4
python cli.py admet --mw 300 --logp 2.5 --hbd 2 --hba 4 --psa 80
python cli.py affinity --delta-g -7 --temperature 298.15
```

`sample_docking.json` demonstrates the expected scoring input. The distinct `sample_payload.json` is an **audit task**, not valid docking data. Ligand properties are optional: `rotatable_bonds`, `mw`, `logp`, `hbd`, `hba`, `psa`. Coordinates are in Å, score in arbitrary kcal/mol-like units, temperature in Kelvin.

## API and audit service

Start the local API:

```bash
python -m uvicorn vina_docking_scorer.server:create_app --factory --host 127.0.0.1 --port 8000
```

Open `http://127.0.0.1:8000/` for the browser calculator or `/docs` for interactive API documentation.

| Endpoint | Purpose |
| --- | --- |
| `GET /health` | Server health |
| `POST /api/score` | Score ligand/receptor JSON (same schema as the example, with optional top-level `temperature`) |
| `POST /api/audit` | Generic operational-threshold audit, separate from docking |
| `POST /api/chat` | Deterministic mock supervisory response |
| `GET /api/audit/logs` | In-memory audit trail and integrity status |
| `GET /metrics` | Simple operational counters (JSON) |

The audit subsystem implements simple rules and a chained HMAC-SHA256 signature log. **It is not an authenticated, persistent, regulatory-grade audit service**: logs and the randomly generated fallback key are lost on restart, and the PHI regex checks are incomplete safeguards rather than guaranteed de-identification or HIPAA compliance. The mock conversational response is not a validated inference service.

For Docker Compose, set a strong random `AUDIT_SECRET_KEY` before starting:

```bash
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_hex(32))')"
docker compose up --build
```

The Compose port is bound to localhost by default. Do not expose the unauthenticated audit API to the public internet. Keep sensitive patient information out of demo requests.

## Development and tests

```bash
python -m pytest -q
node --check web/calculator.mjs
node --test web/test_calculator.mjs
```

GitHub Actions runs Python tests across supported interpreter versions, package-install checks, and Node browser-calculator tests. The [Pages workflow](.github/workflows/pages.yml) publishes the `web/` folder when GitHub Pages is enabled for GitHub Actions. A normal production FastAPI deployment requires a server or container and cannot run on GitHub Pages.

## Privacy and license

The static calculator processes user-supplied JSON entirely in the browser and does not send it to a server. The optional backend audit explicitly sends entered fields to the configured FastAPI service. No external inference provider is used by the current mock implementation.

Licensed under the [MIT License](LICENSE).
