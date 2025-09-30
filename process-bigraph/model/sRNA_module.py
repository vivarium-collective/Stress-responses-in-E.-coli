from process_bigraph.composite import Process
from process_bigraph import register_types, ProcessTypes
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------- Process ----------------------------
class SRNARegulator(Process):
    """
    Minimal sRNA regulation of a target mRNA with stress control.
    State (single 'cell' store): s, m, c, P, S
    """
    config_schema = {
        # 'activator' ~ ArcZ/DsrA/RprA on rpoS; 'repressor' ~ OxyS/CyaR
        'mode': {'_type': 'string', '_default': 'activator'},

        # Hill transcription vs stress S in [0,1]
        'k_tx_s': {'_type': 'map', '_default': {'k0': 0.02, 'kmax': 0.30, 'K': 0.30, 'n': 2.0}},
        'k_tx_m': {'_type': 'map', '_default': {'k0': 0.05, 'kmax': 0.50, 'K': 0.40, 'n': 2.0}},

        # Hfq / pairing kinetics
        'H':       {'_type': 'float', '_default': 1.0},
        'kon':     {'_type': 'float', '_default': 3e-4},
        'koff':    {'_type': 'float', '_default': 0.10},
        'kcleave': {'_type': 'float', '_default': 0.10},
        'krescue': {'_type': 'float', '_default': 0.05},

        # decay
        'kdeg_s':  {'_type': 'float', '_default': 0.20},
        'kdeg_m':  {'_type': 'float', '_default': 0.10},

        # translation control
        # activator: ktl = ktl0 + gamma*c/(Kg+c)
        # repressor: ktl = ktl0 / (1 + eta*c/(Ke+c))
        'ktl0':    {'_type': 'float', '_default': 0.50},
        'gamma':   {'_type': 'float', '_default': 3.0},
        'Kg':      {'_type': 'float', '_default': 50.0},
        'eta':     {'_type': 'float', '_default': 4.0},
        'Ke':      {'_type': 'float', '_default': 50.0},

        # stress-dependent protein stabilization
        'Pdeg0':      {'_type': 'float', '_default': 0.02},
        'sigma_stab': {'_type': 'float', '_default': 0.6},
        'Kp':         {'_type': 'float', '_default': 0.3},
        'np':         {'_type': 'float', '_default': 2.0},

        # optional Rho termination relief by complex
        'rho_p0':   {'_type': 'float', '_default': 0.40},
        'rho_delta':{'_type': 'float', '_default': 0.50},
        'rho_Kc':   {'_type': 'float', '_default': 50.0},
    }

    # one-port API (everything in one map store)
    def inputs(self):  return {'cell': 'map[float]'}
    def outputs(self): return {'cell': 'map[float]'}
    def initial_state(self):
        return {'cell': {'s': 0.0, 'm': 5.0, 'c': 0.0, 'P': 0.0, 'S': 0.0}}

    @staticmethod
    def hill(k0, kmax, K, n, S):
        return k0 + kmax * (S**n) / (K**n + S**n + 1e-12)

    def update(self, state, interval):
        x = state['cell']
        S = float(x.get('S', 0.0))
        cfg = self.config

        # transcription
        k_tx_s = self.hill(**cfg['k_tx_s'], S=S)
        k_tx_m_raw = self.hill(**cfg['k_tx_m'], S=S)

        # Rho termination relief by c
        c = x.get('c', 0.0)
        rho_relief = cfg['rho_delta'] * (c / (cfg['rho_Kc'] + c + 1e-12))
        p_rho = cfg['rho_p0'] * (1.0 - rho_relief)
        k_tx_m = k_tx_m_raw * (1.0 - p_rho)

        # pairing
        s, m = x.get('s', 0.0), x.get('m', 0.0)
        H = cfg['H']
        kon, koff = cfg['kon'], cfg['koff']
        kcleave_eff = cfg['kcleave']
        if cfg['mode'] == 'activator':
            ktl = cfg['ktl0'] + cfg['gamma'] * c / (cfg['Kg'] + c + 1e-12)
        else:  # repressor
            ktl = cfg['ktl0'] / (1.0 + cfg['eta'] * c / (cfg['Ke'] + c + 1e-12))
            kcleave_eff *= 1.5

        # stress-dependent protein degradation
        P = x.get('P', 0.0)
        Pdeg = cfg['Pdeg0'] * (1.0 - cfg['sigma_stab'] * (S**cfg['np']) / (cfg['Kp']**cfg['np'] + S**cfg['np'] + 1e-12))

        # ODEs
        ds = k_tx_s - (kon*H*s*m) + koff*c - cfg['kdeg_s']*s
        dm = k_tx_m - (kon*H*s*m) + koff*c - cfg['kdeg_m']*m
        dc = (kon*H*s*m) - (koff + kcleave_eff + cfg['krescue'])*c
        dP = (ktl * m) - Pdeg*P

        # deltas (Euler)
        return {'cell': {
            's': ds*interval,
            'm': dm*interval,
            'c': dc*interval,
            'P': dP*interval,
            'S': 0.0  # left for driver to overwrite each step
        }}
