import tellurium as te
import matplotlib.pyplot as plt
import numpy as np

# ==============================
# Antimony model with events
# ==============================
ant = r"""
model RpoS_sRNA_core()

  //----- Compartment
  compartment cell = 1;

  //----- Species
  species rpoS_mRNA, RpoS, DsrA, RprA, ArcZ, OxyS;
  species C_DsrA, C_RprA, C_ArcZ, C_OxyS;

  //----- Parameters
  k_tx_rpoS   = 5;
  k_tx_DsrA   = 1;     
  k_tx_RprA   = 1;     
  k_tx_ArcZ   = 1;     
  k_tx_OxyS   = 0.5;

  stress_cold = 1;       
  stress_env  = 1;       
  stress_redx = 1;       
  stress_ox   = 1;       

  d_mRNA  = 0.02;        
  d_RpoS  = 0.001;       
  d_sRNA  = 0.05;        

  d_C_act = 0.01;        
  d_C_oxy = 0.08;        

  k_tl_base = 0.8;       
  k_tl_act  = 3.0;       

  kon_D = 0.005; koff_D = 0.02;     
  kon_R = 0.004; koff_R = 0.02;     
  kon_A = 0.004; koff_A = 0.02;     
  kon_O = 0.006; koff_O = 0.03;     

  //=============================
  // Reactions
  //-----------------------------
  J_tx_rpoS:  -> rpoS_mRNA;  k_tx_rpoS;
  J_tx_DsrA:  -> DsrA;       k_tx_DsrA*stress_cold;
  J_tx_RprA:  -> RprA;       k_tx_RprA*stress_env;
  J_tx_ArcZ:  -> ArcZ;       k_tx_ArcZ*stress_redx;
  J_tx_OxyS:  -> OxyS;       k_tx_OxyS*stress_ox;

  J_deg_mRNA: rpoS_mRNA -> ; d_mRNA * rpoS_mRNA;
  J_deg_DsrA: DsrA      -> ; d_sRNA * DsrA;
  J_deg_RprA: RprA      -> ; d_sRNA * RprA;
  J_deg_ArcZ: ArcZ      -> ; d_sRNA * ArcZ;
  J_deg_OxyS: OxyS      -> ; d_sRNA * OxyS;
  J_deg_RpoS: RpoS      -> ; d_RpoS * RpoS;

  J_bind_D: rpoS_mRNA + DsrA -> C_DsrA; kon_D*rpoS_mRNA*DsrA;
  J_unbind_D: C_DsrA -> rpoS_mRNA + DsrA; koff_D*C_DsrA;

  J_bind_R: rpoS_mRNA + RprA -> C_RprA; kon_R*rpoS_mRNA*RprA;
  J_unbind_R: C_RprA -> rpoS_mRNA + RprA; koff_R*C_RprA;

  J_bind_A: rpoS_mRNA + ArcZ -> C_ArcZ; kon_A*rpoS_mRNA*ArcZ;
  J_unbind_A: C_ArcZ -> rpoS_mRNA + ArcZ; koff_A*C_ArcZ;

  J_bind_O: rpoS_mRNA + OxyS -> C_OxyS; kon_O*rpoS_mRNA*OxyS;
  J_unbind_O: C_OxyS -> rpoS_mRNA + OxyS; koff_O*C_OxyS;

  J_deg_Cact: C_DsrA -> ; d_C_act*C_DsrA;
  J_deg_Crpr: C_RprA -> ; d_C_act*C_RprA;
  J_deg_Carc: C_ArcZ -> ; d_C_act*C_ArcZ;
  J_deg_Coxy: C_OxyS -> ; d_C_oxy*C_OxyS;

  J_tl_free:  -> RpoS; k_tl_base * rpoS_mRNA;
  J_tl_D:     -> RpoS; k_tl_act  * C_DsrA;
  J_tl_R:     -> RpoS; k_tl_act  * C_RprA;
  J_tl_A:     -> RpoS; k_tl_act  * C_ArcZ;

  //=============================
  // Initial conditions
  rpoS_mRNA = 0;
  RpoS      = 0;
  DsrA      = 0;  
  RprA      = 0;  
  ArcZ      = 0;  
  OxyS      = 0;
  C_DsrA    = 0;  
  C_RprA    = 0;  
  C_ArcZ    = 0;  
  C_OxyS    = 0;

  //=============================
  // Events: oxidative stress pulse
  at (time > 200): stress_ox = 4;
  at (time > 400): stress_ox = 1;

end
"""

# ==============================
# Load model
# ==============================
r = te.loada(ant)
species = r.getFloatingSpeciesIds()
RpoS_idx = species.index('RpoS') + 1

# ==============================
# Run simulation with pulse
# ==============================
pulse = r.simulate(0, 600, 601)
t = pulse[:, 0]

plt.figure(figsize=(8, 5))
plt.plot(t, pulse[:, RpoS_idx], label="RpoS (protein)")
plt.axvspan(200, 400, color="red", alpha=0.2, label="Oxidative stress pulse")
plt.xlabel("Time");
plt.ylabel("RpoS (a.u.)")
plt.title("RpoS dynamics under oxidative stress pulse (200–400)")
plt.legend();
plt.tight_layout();
# plt.show()
print(pulse)