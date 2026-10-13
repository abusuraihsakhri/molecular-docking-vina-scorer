/**
 * Browser implementation of the repository's simplified Vina-like scoring
 * and descriptor heuristics. No molecular docking search is performed.
 * No data leaves the browser.
 */
const weights = Object.freeze({
  gauss1: -0.0356, gauss2: -0.00516, repulsion: 0.840,
  hydrophobic: -0.0351, hbond: -0.587, torsional: 0.0585
});
const round4 = x => Math.round((x + Number.EPSILON) * 1e4) / 1e4;
function finite(value, name) {
  if (typeof value !== "number" || !Number.isFinite(value)) {
    throw new Error(`${name} must be a finite number`);
  }
  return value;
}
function nonNegativeInt(value, name) {
  if (!Number.isInteger(value) || value < 0) {
    throw new Error(`${name} must be a non-negative integer`);
  }
  return value;
}
function validatedAtoms(atoms, name) {
  if (!Array.isArray(atoms) || atoms.length < 1 || atoms.length > 128) {
    throw new Error(`${name} must contain 1–128 atoms`);
  }
  for (const [index, a] of atoms.entries()) {
    if (!a || typeof a.element !== "string" || !a.element.trim()) {
      throw new Error(`${name}[${index}] needs an element`);
    }
    for (const key of ["x", "y", "z"]) finite(a[key], `${name}[${index}].${key}`);
  }
  return atoms;
}
export function score(ligandAtoms, receptorAtoms, rotatable = 0) {
  validatedAtoms(ligandAtoms, "ligand.atoms");
  validatedAtoms(receptorAtoms, "receptor.atoms");
  nonNegativeInt(rotatable, "rotatable_bonds");
  const raw = { gauss1: 0, gauss2: 0, repulsion: 0, hydrophobic: 0, hbond: 0 };
  for (const a of ligandAtoms) {
    for (const b of receptorAtoms) {
      const d = Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
      raw.gauss1 += Math.exp(-((d - 3) / .5) ** 2);
      raw.gauss2 += Math.exp(-((d - 3) / 2) ** 2);
      if (d < 2) raw.repulsion += (2 - d) ** 2;
      if (a.hydrophobic && b.hydrophobic) {
        raw.hydrophobic += d <= 1.5 ? 1 : d >= 4 ? 0 : 1 - (d - 1.5) / 2.5;
      }
      if ((a.donor && b.acceptor) || (a.acceptor && b.donor)) {
        raw.hbond += Math.exp(-((d - 2.8) / .7) ** 2);
      }
    }
  }
  const total = raw.gauss1 * weights.gauss1 + raw.gauss2 * weights.gauss2 +
    raw.repulsion * weights.repulsion + raw.hydrophobic * weights.hydrophobic +
    raw.hbond * weights.hbond + rotatable * weights.torsional;
  return {
    gauss1: round4(raw.gauss1 * weights.gauss1),
    gauss2: round4(raw.gauss2 * weights.gauss2),
    repulsion: round4(raw.repulsion * weights.repulsion),
    hydrophobic: round4(raw.hydrophobic * weights.hydrophobic),
    hbond: round4(raw.hbond * weights.hbond),
    torsional: round4(rotatable * weights.torsional),
    total: round4(total)
  };
}
export function affinity(deltaG, temperature = 298.15) {
  finite(deltaG, "delta_g");
  finite(temperature, "temperature");
  if (temperature <= 0) throw new Error("temperature must be positive");
  const exponent = deltaG / (0.001987 * temperature);
  if (exponent < -700 || exponent > 680) throw new Error("Energy/temperature combination is outside the numeric range");
  const ki = Math.exp(exponent);
  const ki_nm = ki * 1e9;
  const classification = ki_nm < 10 ? "very_strong" : ki_nm < 100 ? "strong" :
    ki_nm < 1000 ? "moderate" : ki_nm < 10000 ? "weak" : "very_weak";
  return { delta_g: round4(deltaG), ki, ki_nm: round4(ki_nm), kd: ki, classification };
}
export function lipinski(mw, logp, hbd, hba) {
  finite(mw, "mw"); finite(logp, "logp");
  nonNegativeInt(hbd, "hbd"); nonNegativeInt(hba, "hba");
  if (mw <= 0) throw new Error("mw must be positive");
  const details = [];
  if (mw > 500) details.push(`MW ${mw.toFixed(1)} > 500`);
  if (logp > 5) details.push(`LogP ${logp.toFixed(2)} > 5`);
  if (hbd > 5) details.push(`HBD ${hbd} > 5`);
  if (hba > 10) details.push(`HBA ${hba} > 10`);
  return { mw, logp, hbd, hba, passes: !details.length, violations: details.length, violation_details: details };
}
export function admet(mw, logp, hbd, hba, psa = 80, rotatable = 3) {
  finite(mw,"mw"); finite(logp,"logp"); finite(psa,"psa");
  nonNegativeInt(hbd,"hbd"); nonNegativeInt(hba,"hba"); nonNegativeInt(rotatable,"rotatable");
  if (mw <= 0 || psa < 0) throw new Error("mw must be positive and psa non-negative");
  const absorption = logp < 0 || logp > 5 || psa > 140 || mw > 500 ? "low" :
    logp >= 1 && logp <= 3 && psa < 100 ? "high" : "moderate";
  const distribution = logp < 0 ? "low" : logp > 3 ? "high" : "moderate";
  const metabolism = mw > 500 || hba > 8 ? "high" : mw < 300 && hba < 4 ? "low" : "moderate";
  const excretion = mw > 500 || psa > 120 ? "low" : mw < 300 && psa < 80 ? "high" : "moderate";
  const risk = Number(logp > 5) + Number(mw > 500) + Number(hba > 10) + Number(psa < 20);
  const toxicity_risk = risk >= 3 ? "high" : risk >= 1 ? "moderate" : "low";
  const bioScore = Number(mw >= 200 && mw <= 500) + Number(logp >= 0 && logp <= 3) +
    Number(hbd <= 5) + Number(hba <= 10) + Number(psa <= 140) + Number(rotatable <= 10);
  const oral_bioavailability = bioScore >= 5 ? "high" : bioScore >= 3 ? "moderate" :
    bioScore >= 1 ? "low" : "very_low";
  return { absorption, distribution, metabolism, excretion, toxicity_risk,
    oral_bioavailability, details: { mw, logp, hbd, hba, psa, num_rotatable: rotatable, bio_score: bioScore } };
}
export function analyze(input) {
  if (!input || !input.ligand || !input.receptor) {
    throw new Error("JSON must include ligand and receptor objects");
  }
  const p = input.ligand.properties ?? {};
  const rotatable = input.rotatable_bonds ?? p.rotatable_bonds ?? 0;
  const energy = score(input.ligand.atoms, input.receptor.atoms, rotatable);
  const out = { score: energy, binding_affinity: affinity(energy.total, input.temperature ?? 298.15), lipinski: null, admet: null };
  if (p.mw != null) {
    out.lipinski = lipinski(p.mw, p.logp ?? 0, p.hbd ?? 0, p.hba ?? 0);
    out.admet = admet(p.mw, p.logp ?? 0, p.hbd ?? 0, p.hba ?? 0, p.psa ?? 80, rotatable);
  }
  return out;
}
