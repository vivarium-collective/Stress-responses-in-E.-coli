"""
 Tellurium Antimony model for sRNA regulation of RpoS


Core biology:
- RpoS mRNA and protein
- sRNAs: DsrA, RprA, ArcZ (activators), OxyS (repressor)
- sRNA–mRNA complexes
- RNAP–sigma competition: E + σ70 <-> E70 ; E + σS <-> ES
- RpoS transcription depends on RNAP·σS (ES)
"""

import tellurium as te
import matplotlib.pyplot as plt

# -----------------------------
# Antimony model definition
# -----------------------------
antimony_str = r"""
model RpoS_sRNA_competition()

  // ===== Species
  species rpoS_mRNA, RpoS, DsrA, RprA, ArcZ, OxyS;
  species C_DsrA, C_RprA, C_ArcZ, C_OxyS;
  species E70, ES;   // RNAP holoenzymes

  // ===== Parameters
  E_tot = 400; Sig70_tot = 800; SigS_tot = 150;
  kon70 = 0.002; koff70 = 0.02;
  konS  = 0.003; koffS  = 0.02;

  k_tx_rpoS_max = 6;  K_RNAP = 40;

  k_tx_DsrA = 1; k_tx_RprA = 1; k_tx_ArcZ = 1; k_tx_OxyS = 0.5;

  stress_cold = 1; stress_env = 1; stress_redx = 1; stress_ox = 1;

  d_mRNA=0.02; d_RpoS=0.001; d_sRNA=0.05;
  d_C_act=0.01; d_C_oxy=0.08;

  k_tl_base=0.8; k_tl_act=3.0;

  kon_D=0.005; koff_D=0.02;
  kon_R=0.004; koff_R=0.02;
  kon_A=0.004; koff_A=0.02;
  kon_O=0.006; koff_O=0.03;

  // ===== Assignment rules for free pools
  E_free     := E_tot - E70 - ES;
  Sig70_free := Sig70_tot - E70;
  SigS_free  := SigS_tot - ES;

  // ===== RNAP binding
  R_bind70: E_free + Sig70_free -> E70; kon70*E_free*Sig70_free;
  R_unbd70: E70 -> E_free + Sig70_free; koff70*E70;

  R_bindS: E_free + SigS_free -> ES; konS*E_free*SigS_free;
  R_unbdS: ES -> E_free + SigS_free; koffS*ES;

  // ===== Transcription
  J_tx_rpoS: -> rpoS_mRNA; k_tx_rpoS_max * ES/(ES + K_RNAP);
  J_tx_DsrA: -> DsrA; k_tx_DsrA*stress_cold;
  J_tx_RprA: -> RprA; k_tx_RprA*stress_env;
  J_tx_ArcZ: -> ArcZ; k_tx_ArcZ*stress_redx;
  J_tx_OxyS: -> OxyS; k_tx_OxyS*stress_ox;

  // ===== Decay
  J_deg_mRNA: rpoS_mRNA -> ; d_mRNA*rpoS_mRNA;
  J_deg_DsrA: DsrA -> ; d_sRNA*DsrA;
  J_deg_RprA: RprA -> ; d_sRNA*RprA;
  J_deg_ArcZ: ArcZ -> ; d_sRNA*ArcZ;
  J_deg_OxyS: OxyS -> ; d_sRNA*OxyS;
  J_deg_RpoS: RpoS -> ; d_RpoS*RpoS;

  // ===== sRNA binding/unbinding
  J_bind_D: rpoS_mRNA + DsrA -> C_DsrA; kon_D*rpoS_mRNA*DsrA;
  J_unbd_D: C_DsrA -> rpoS_mRNA + DsrA; koff_D*C_DsrA;

  J_bind_R: rpoS_mRNA + RprA -> C_RprA; kon_R*rpoS_mRNA*RprA;
  J_unbd_R: C_RprA -> rpoS_mRNA + RprA; koff_R*C_RprA;

  J_bind_A: rpoS_mRNA + ArcZ -> C_ArcZ; kon_A*rpoS_mRNA*ArcZ;
  J_unbd_A: C_ArcZ -> rpoS_mRNA + ArcZ; koff_A*C_ArcZ;

  J_bind_O: rpoS_mRNA + OxyS -> C_OxyS; kon_O*rpoS_mRNA*OxyS;
  J_unbd_O: C_OxyS -> rpoS_mRNA + OxyS; koff_O*C_OxyS;

  // ===== Complex decay
  J_deg_Cact: C_DsrA -> ; d_C_act*C_DsrA;
  J_deg_Crpr: C_RprA -> ; d_C_act*C_RprA;
  J_deg_Carc: C_ArcZ -> ; d_C_act*C_ArcZ;
  J_deg_Coxy: C_OxyS -> ; d_C_oxy*C_OxyS;

  // ===== Translation
  J_tl_free: -> RpoS; k_tl_base*rpoS_mRNA;
  J_tl_D:    -> RpoS; k_tl_act*C_DsrA;
  J_tl_R:    -> RpoS; k_tl_act*C_RprA;
  J_tl_A:    -> RpoS; k_tl_act*C_ArcZ;

  // ===== Initial conditions
  rpoS_mRNA=0; RpoS=0; DsrA=0; RprA=0; ArcZ=0; OxyS=0;
  C_DsrA=0; C_RprA=0; C_ArcZ=0; C_OxyS=0;
  E70=0; ES=0;

end
"""

# -----------------------------
# Helper functions
# -----------------------------
def load_model():
    return te.loada(antimony_str)

def simulate_baseline(t_end=600, points=601):
    r = load_model()
    return r.simulate(0, t_end, points)

def simulate_stress(stress_param, value, t_end=600, points=601):
    r = load_model()
    r[stress_param] = value
    return r.simulate(0, t_end, points)

def plot_species(result, species_list, title="Simulation"):
    t = result[:,0]
    plt.figure(figsize=(7,4))
    for sid in species_list:
        idx = r.getFloatingSpeciesIds().index(sid) + 1
        plt.plot(t, result[:,idx], label=sid)
    plt.xlabel("Time"); plt.ylabel("Concentration (a.u.)")
    plt.title(title); plt.legend(); plt.tight_layout()
    plt.show()

# -----------------------------
# Example run
# -----------------------------
if __name__ == "__main__":
    r = load_model()
    res = r.simulate(0, 600, 601)
    r.plot(res)  # quick plot
