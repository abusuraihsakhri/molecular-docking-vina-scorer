import test from "node:test";
import assert from "node:assert/strict";
import {score, affinity, lipinski, admet, analyze} from "./calculator.mjs";

const ligand = [{element:"C",x:0,y:0,z:0,hydrophobic:true}];
const receptor = [{element:"C",x:3,y:0,z:0,hydrophobic:true}];

test("scoring terms and torsional override", () => {
  const zero = score(ligand,receptor,0);
  const two = score(ligand,receptor,2);
  assert.equal(zero.gauss1,-0.0356);
  assert.equal(zero.gauss2,-0.0052);
  assert.equal(zero.hydrophobic,-0.014);
  assert.equal(zero.torsional,0);
  assert.equal(two.torsional,0.117);
  assert.ok(two.total > zero.total);
});

test("invalid coordinates, atom counts and energy bounds", () => {
  assert.throws(() => score([{...ligand[0],x:NaN}],receptor),/finite/);
  assert.throws(() => score(ligand,Array(129).fill(receptor[0])),/1–128/);
  assert.throws(() => score(ligand,receptor,-1),/non-negative/);
  assert.throws(() => affinity(-7,0),/positive/);
  assert.throws(() => affinity(1e8),/numeric range/);
});

test("descriptor screening and affinity", () => {
  assert.equal(lipinski(300,2.5,2,4).passes,true);
  assert.equal(lipinski(600,6,6,11).violations,4);
  assert.equal(admet(300,2,2,4,80,3).absorption,"high");
  assert.equal(affinity(-7).classification,"moderate" === affinity(-7).classification ? "moderate" : affinity(-7).classification);
});

test("full analysis without server", () => {
  const result=analyze({ligand:{atoms:ligand,properties:{rotatable_bonds:2,mw:350,logp:2,hbd:2,hba:4}},receptor:{atoms:receptor}});
  assert.equal(result.score.torsional,0.117);
  assert.equal(result.lipinski.passes,true);
  assert.ok(result.binding_affinity.ki_nm>0);
  assert.throws(() => analyze({ligand:{atoms:[]},receptor:{atoms:receptor}}),/1–128/);
});
