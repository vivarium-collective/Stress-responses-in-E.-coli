import tellurium as te
import matplotlib.pyplot as plt
import numpy as np

# =========================================
# Antimony model with RNAP–sigma competition
# =========================================
ant = r"""
model RpoS_sRNA_competition()

  //=============================
  // Compartment
  compartment cell = 1;

  //=============================
  // Core RpoS/sRNA species (same as before)
  species rpoS_mRNA, RpoS, DsrA, RprA, ArcZ, OxyS;
  species C_DsrA, C_RprA, C_ArcZ, C_OxyS;

  //=============================
  // RNAP–sigma competition block
  // Species for holoenzymes (tracked explicitly)
  species E70, ES;        // RNAP·σ70 and RNAP·σS holoenzymes

  // Totals (parameters)
  E_tot       = 400.0;    // total RNAP core
  Sig70_tot   = 800.0;    // total σ70 pool
  SigS_tot    = 150.0;    // total σS pool (RpoS)
  // Binding kinetics
  kon70 = 0.002;  koff70 = 0.02;   // E + σ70 <-> E70
  konS  = 0.003;  koffS  = 0.02;   // E + σS  <-> ES (slightly stronger here)

  // Assignment (mass balance) for *free* species
  // These are *not* independent species; they’re computed each step.
  E_free      := E_tot    - E70 - ES;
  Sig70_free  := Sig70_tot - E70;
  SigS_free   := SigS_tot  - ES;

  // Reversible binding to form holoenzymes
  R_bind70: E_free + Sig70_free -> E70;  kon70 * E_free * Sig70_free;
  R_unbd70: E70 -> E_free + Sig70_free;  koff70 * E70;

  R_bindS:  E_free + SigS_free  -> ES;   konS  * E_free * SigS_free;
  R_unbdS:  ES   -> E_free + SigS_free;  koffS  * ES;

  //=============================
  // Transcription capacities
  // RpoS transcription is limited by ES availability (saturable)
  k_tx_rpoS_max = 6.0;   // max rate if ES is abundant
  K_RNAP        = 40.0;  // ES level for half-max RpoS transcription

  // sRNA baseline transcription (before stress scaling)
  k_tx_DsrA   = 1;
  k_tx_RprA   = 1;
  k_tx_ArcZ   = 1;
  k_tx_OxyS   = 0.5;

  // Stress multipliers (editable at runtime/events)
  stress_cold = 1;       // ↑DsrA
  stress_env  = 1;       // ↑RprA
  stress_redx = 1;       // ↑ArcZ
  stress_ox   = 1;       // ↑OxyS

  // Decay rates
  d_mRNA  = 0.02;
  d_RpoS  = 0.001;
  d_sRNA  = 0.05;

  // Complex-specific decay
  d_C_act = 0.01;
  d_C_oxy = 0.08;

  // Translation
  k_tl_base = 0.8;
  k_tl_act  = 3.0;

  // sRNA–mRNA binding/unbinding
  kon_D = 0.005; koff_D = 0.02;
  kon_R = 0.004; koff_R = 0.02;
  kon_A = 0.004; koff_A = 0.02;
  kon_O = 0.006; koff_O = 0.03;

  //=============================
  // Reactions
  //-----------------------------
  // NOTE: RpoS transcription now depends on ES (RNAP·σS)
  J_tx_rpoS:  -> rpoS_mRNA;  k_tx_rpoS_max * ES / (ES + K_RNAP);

  // sRNA transcription (stress-scaled, keep as before)
  J_tx_DsrA:  -> DsrA;       k_tx_DsrA*stress_cold;
  J_tx_RprA:  -> RprA;       k_tx_RprA*stress_env;
  J_tx_ArcZ:  -> ArcZ;       k_tx_ArcZ*stress_redx;
  J_tx_OxyS:  -> OxyS;       k_tx_OxyS*stress_ox;

  // Decay
  J_deg_mRNA: rpoS_mRNA -> ; d_mRNA * rpoS_mRNA;
  J_deg_DsrA: DsrA      -> ; d_sRNA * DsrA;
  J_deg_RprA: RprA      -> ; d_sRNA * RprA;
  J_deg_ArcZ: ArcZ      -> ; d_sRNA * ArcZ;
  J_deg_OxyS: OxyS      -> ; d_sRNA * OxyS;
  J_deg_RpoS: RpoS      -> ; d_RpoS * RpoS;

  // Binding/unbinding (activators)
  J_bind_D: rpoS_mRNA + DsrA -> C_DsrA;           kon_D*rpoS_mRNA*DsrA;
  J_unbd_D: C_DsrA -> rpoS_mRNA + DsrA;           koff_D*C_DsrA;

  J_bind_R: rpoS_mRNA + RprA -> C_RprA;           kon_R*rpoS_mRNA*RprA;
  J_unbd_R: C_RprA -> rpoS_mRNA + RprA;           koff_R*C_RprA;

  J_bind_A: rpoS_mRNA + ArcZ -> C_ArcZ;           kon_A*rpoS_mRNA*ArcZ;
  J_unbd_A: C_ArcZ -> rpoS_mRNA + ArcZ;           koff_A*C_ArcZ;

  // Binding/unbinding (repressor)
  J_bind_O: rpoS_mRNA + OxyS -> C_OxyS;           kon_O*rpoS_mRNA*OxyS;
  J_unbd_O: C_OxyS -> rpoS_mRNA + OxyS;           koff_O*C_OxyS;

  // Complex decay
  J_deg_Cact: C_DsrA -> ; d_C_act*C_DsrA;
  J_deg_Crpr: C_RprA -> ; d_C_act*C_RprA;
  J_deg_Carc: C_ArcZ -> ; d_C_act*C_ArcZ;
  J_deg_Coxy: C_OxyS -> ; d_C_oxy*C_OxyS;

  // Translation
  J_tl_free:  -> RpoS; k_tl_base * rpoS_mRNA;
  J_tl_D:     -> RpoS; k_tl_act  * C_DsrA;
  J_tl_R:     -> RpoS; k_tl_act  * C_RprA;
  J_tl_A:     -> RpoS; k_tl_act  * C_ArcZ;

  //=============================
  // Initial conditions
  rpoS_mRNA = 0;
  RpoS      = 0;
  DsrA      = 0;  RprA = 0;  ArcZ = 0;  OxyS = 0;
  C_DsrA    = 0;  C_RprA = 0; C_ArcZ = 0; C_OxyS = 0;

  // Start with no pre-formed holoenzymes
  E70 = 0;
  ES  = 0;

end
"""

# Load the model
r = te.loada(ant)
species = r.getFloatingSpeciesIds()
RpoS_idx = species.index('RpoS') + 1

# --------------------------
# 1) Baseline with competition
# --------------------------
base = r.simulate(0, 800, 801)
t = base[:, 0]

plt.figure(figsize=(7.5,4.2))
plt.plot(t, base[:, RpoS_idx])
plt.xlabel("Time"); plt.ylabel("RpoS (a.u.)")
plt.title("Baseline with RNAP–σ competition (no stress)")
plt.tight_layout(); plt.show()

# --------------------------
# 2) Show effect of σ70 competition vs σS on RpoS output
#    Sweep Sig70_tot upward: more σ70 steals RNAP from σS → less ES → lower RpoS tx
# --------------------------
def steady_RpoS_for_sigma70(s70_total):
    r.resetAll()
    r["Sig70_tot"] = float(s70_total)
    out = r.simulate(0, 1200, 1201)
    return out[-1, RpoS_idx]

sigma70_scan = np.linspace(200, 1400, 13)   # try a wide range
RpoS_steady  = [steady_RpoS_for_sigma70(x) for x in sigma70_scan]

plt.figure(figsize=(7.5,4.2))
plt.plot(sigma70_scan, RpoS_steady, marker="o")
plt.xlabel("σ70 total (molecules)")
plt.ylabel("Steady RpoS (a.u.)")
plt.title("Sigma-factor competition curve: increasing σ70 suppresses RpoS")
plt.tight_layout(); plt.show()

# --------------------------
# 3)  two time-courses at different σ70_tot
# --------------------------
def timecourse_with_sigma70(s70_total):
    r.resetAll()
    r["Sig70_tot"] = float(s70_total)
    return r.simulate(0, 800, 801)

tc_low  = timecourse_with_sigma70(400)   # less σ70 → more ES → stronger RpoS tx
tc_high = timecourse_with_sigma70(1200)  # more σ70 → less ES → weaker RpoS tx

plt.figure(figsize=(7.5,4.2))
plt.plot(tc_low[:,0],  tc_low[:, RpoS_idx],  label="σ70_tot = 400")
plt.plot(tc_high[:,0], tc_high[:, RpoS_idx], label="σ70_tot = 1200")
plt.xlabel("Time"); plt.ylabel("RpoS (a.u.)")
plt.title("Time-courses under different σ70 pools")
plt.legend()
plt.tight_layout(); plt.show()
