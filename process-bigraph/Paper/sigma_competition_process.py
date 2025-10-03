
"""
SigmaCompetition Process and Composite utilities (Mauri & Klumpp, 2014 style).

This module:
- Defines a process-bigraph Process: SigmaCompetition.
- Provides build_core, build_alloc_composite, and step_alloc_once helpers.
- Robustly normalizes Composite.update results (dict or list-of-dicts).
"""

from __future__ import annotations
from typing import Dict, Iterable, Mapping, MutableMapping, Tuple

import numpy as np
from scipy.optimize import fsolve

from process_bigraph import register_types, ProcessTypes
from process_bigraph.composite import Process, Composite


# =============================================================================
# Process: SigmaCompetition
# =============================================================================
class SigmaCompetition(Process):
    """
    Steady-state RNAP allocation among sigma factors (E_free, E·σ70, E·σAlt),
    plus simple promoter-like outputs (J_σ70, J_σAlt).
    """

    config_schema = {
        # Totals (absolute or a.u.; panel C may use μM-like units consistently)
        'RNAP_total': {'_type': 'float', '_default': 11400.0},
        'sigma70_total': {'_type': 'float', '_default': 5700.0},
        'sigmaS_total': {'_type': 'float', '_default': 2000.0},

        # Effective dissociation constants (lower => tighter binding)
        'Kd_sigma70': {'_type': 'float', '_default': 1.0},
        'Kd_sigmaS':  {'_type': 'float', '_default': 20.0},

        # Simple promoter model J = n * a * E_sigma / (K + E_sigma)
        'K_prom': {'_type': 'float', '_default': 100.0},
        'a_prom': {'_type': 'float', '_default': 1.0},
        'n_promoters_sigma70': {'_type': 'integer', '_default': 200},
        'n_promoters_sigmaS':  {'_type': 'integer', '_default': 200},
    }

    def inputs(self) -> Mapping[str, str]:
        return {}

    def outputs(self) -> Mapping[str, str]:
        return {
            'E_free': 'float',
            'E_sigma70': 'float',
            'E_sigmaS': 'float',
            'J_sigma70': 'float',
            'J_sigmaS': 'float',
        }

    # ------------------------ Core equations ------------------------
    @staticmethod
    def _equations(
        vars_: Iterable[float],
        RNAP_total: float,
        s70: float,
        sS: float,
        Kd70: float,
        KdS: float,
    ) -> Tuple[float, float, float]:
        E_free, E_sigma70, E_sigmaS = vars_
        eps = 1e-12
        return (
            E_free + E_sigma70 + E_sigmaS - RNAP_total,                 # RNAP conservation
            E_sigma70 - ((s70 - E_sigma70) * E_free / (Kd70 + eps)),    # σ70 binding
            E_sigmaS  - ((sS  - E_sigmaS)  * E_free / (KdS  + eps)),    # alt sigma binding
        )

    def _solve_allocation(
        self,
        RNAP_total: float,
        sigma70_total: float,
        sigmaS_total: float,
        Kd_sigma70: float,
        Kd_sigmaS: float,
    ) -> Tuple[float, float, float]:
        # Reasonable initial guess
        guess = [max(RNAP_total * 0.5, 1.0), RNAP_total * 0.25, RNAP_total * 0.25]
        sol = fsolve(
            self._equations,
            guess,
            args=(RNAP_total, sigma70_total, sigmaS_total, Kd_sigma70, Kd_sigmaS),
            xtol=1e-10,
            maxfev=2000,
        )
        E_free, E70, ES = [max(float(x), 0.0) for x in sol]

        # Clamp to conservation
        total = E_free + E70 + ES
        if total > 0:
            scale = RNAP_total / total
            E_free *= scale
            E70    *= scale
            ES     *= scale
        return E_free, E70, ES

    @staticmethod
    def _promoter_rate(E_sigma: float, K_prom: float, a_prom: float, n_promoters: int) -> float:
        eps = 1e-12
        return float(n_promoters) * float(a_prom) * (E_sigma / (K_prom + E_sigma + eps))

    def update(self, state: Mapping, interval: float) -> Dict[str, float]:
        cfg = self.config
        E_free, E70, ES = self._solve_allocation(
            RNAP_total=float(cfg['RNAP_total']),
            sigma70_total=float(cfg['sigma70_total']),
            sigmaS_total=float(cfg['sigmaS_total']),
            Kd_sigma70=float(cfg['Kd_sigma70']),
            Kd_sigmaS=float(cfg['Kd_sigmaS']),
        )
        J70 = self._promoter_rate(E70, cfg['K_prom'], cfg['a_prom'], cfg['n_promoters_sigma70'])
        JS  = self._promoter_rate(ES,  cfg['K_prom'], cfg['a_prom'], cfg['n_promoters_sigmaS'])
        return {
            'E_free': E_free, 'E_sigma70': E70, 'E_sigmaS': ES,
            'J_sigma70': J70, 'J_sigmaS': JS
        }


# =============================================================================
# Composite helpers
# =============================================================================
_EXPECTED_KEYS = ('E_free', 'E_sigma70', 'E_sigmaS', 'J_sigma70', 'J_sigmaS')


def build_core():
    """Return a fresh core with SigmaCompetition registered."""
    core = register_types(ProcessTypes())
    core.register_process("SigmaCompetition", SigmaCompetition)
    return core


def build_alloc_composite(core, config: Mapping[str, float]) -> Composite:
    """
    Build a Composite with one SigmaCompetition node, wired to top-level stores.
    """
    spec = {
        'E_free':    {'_type': 'float', '_value': 0.0},
        'E_sigma70': {'_type': 'float', '_value': 0.0},
        'E_sigmaS':  {'_type': 'float', '_value': 0.0},
        'J_sigma70': {'_type': 'float', '_value': 0.0},
        'J_sigmaS':  {'_type': 'float', '_value': 0.0},
        'alloc': {
            '_type': 'process',
            'address': 'local:SigmaCompetition',
            'config': dict(config),
            '_outputs': {
                'E_free': 'float',
                'E_sigma70': 'float',
                'E_sigmaS': 'float',
                'J_sigma70': 'float',
                'J_sigmaS': 'float',
            },
            'outputs': {
                'E_free':    ['E_free'],
                'E_sigma70': ['E_sigma70'],
                'E_sigmaS':  ['E_sigmaS'],
                'J_sigma70': ['J_sigma70'],
                'J_sigmaS':  ['J_sigmaS'],
            },
        },
    }
    return Composite(spec, core=core)


def _collect_numbers(obj, out: MutableMapping[str, float]) -> None:
    """Recursively collect numeric leaves matching expected keys."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in _EXPECTED_KEYS and isinstance(v, (int, float)):
                out[k] = out.get(k, 0.0) + float(v)
            else:
                _collect_numbers(v, out)
    elif isinstance(obj, (list, tuple)):
        for item in obj:
            _collect_numbers(item, out)


def normalize_updates(raw) -> Dict[str, float]:
    """Normalize Composite.update return into a flat dict (keys in _EXPECTED_KEYS)."""
    flat: Dict[str, float] = {}
    _collect_numbers(raw, flat)
    return flat


def step_alloc_once(core, config: Mapping[str, float]) -> Dict[str, float]:
    """
    Build Composite, call update once, and return the resulting values, applied
    as deltas to zero-initialized stores. Falls back to direct process.update
    if Composite emits nothing.
    """
    comp = build_alloc_composite(core, config)
    state: Dict[str, float] = {k: 0.0 for k in _EXPECTED_KEYS}

    raw = comp.update(state={}, interval=1.0)
    deltas = normalize_updates(raw)

    if not deltas:
        proc = SigmaCompetition(core=core, config=dict(config))
        direct = proc.update(state={}, interval=1.0)
        deltas = {k: float(direct.get(k, 0.0)) for k in _EXPECTED_KEYS}

    for k, dv in deltas.items():
        state[k] = state.get(k, 0.0) + float(dv)
    for k in _EXPECTED_KEYS:
        state.setdefault(k, 0.0)
    return state