"""Datasets used in the article: wines (sensory), butters and margarines (spectral), and the four
simulated datasets."""
from __future__ import annotations

from importlib import resources

import numpy as np

from .core import prep_blocks, snv

__all__ = ["load_wines", "load_butters", "load_simulated"]


def _read(name, header=True):
    with resources.files("wecops").joinpath("data", name).open("r", encoding="utf-8") as fh:
        lines = fh.read().strip().splitlines()
    if not header:
        return None, None, np.array([[float(v) for v in l.split(",")] for l in lines])
    cols = lines[0].split(",")[1:]
    rows = [l.split(",") for l in lines[1:]]
    return cols, [r[0] for r in rows], np.array([[float(v) for v in r[1:]] for r in rows])


def load_wines(prep=True):
    """Free-choice profiling of 8 wines by 4 judges (4, 3, 4 and 3 personal descriptors scored
    from 0 to 10); Williams & Langron (1984), as reproduced in Hanafi & Qannari (2008).

    Returns a dict with ``blocks`` (list of 4 arrays, centred and of unit inertia if ``prep``),
    ``variables``, ``individuals`` and ``block_names``. The article uses R = 3.
    """
    blocks, variables = [], []
    for k in range(1, 5):
        cols, rows, X = _read(f"wines_judge{k}.csv")
        blocks.append(X)
        variables.append(cols)
    return dict(blocks=prep_blocks(blocks) if prep else blocks, variables=variables, individuals=rows,
                block_names=[f"judge {k}" for k in range(1, 5)])


def load_butters(prep=True):
    """Visible (VIS, 200 variables), near infrared (NIR, 850) and mid infrared (MIR, 1816) spectra
    of 21 fat products, 8 butters and 13 margarines; Jaillais, Morrin & Downey (2007).

    With ``prep``, each spectrum is SNV-transformed, then each block is centred and scaled to unit
    inertia, as in the article (R = 4, and R = 3).
    """
    blocks, variables, rows = [], [], None
    for b in ("VIS", "NIR", "MIR"):
        cols, rows, X = _read(f"butters_{b}.csv")
        blocks.append(snv(X) if prep else X)
        variables.append(cols)
    return dict(blocks=prep_blocks(blocks) if prep else blocks, variables=variables, individuals=rows,
                block_names=["VIS", "NIR", "MIR"])


_SIM = {1: 3, 2: 6, 3: 3, 4: 3}
_SIM_R = {1: 5, 2: 10, 3: 5, 4: 5}


def load_simulated(case):
    """Simulated datasets 1 to 4 of the article (blocks already centred and of unit inertia).

    1: three blocks, one-dimensional structures (R = 5); 2: six blocks of 10 to 40 variables,
    one-dimensional structures (R = 10); 3: one block with surplus proper structures (R = 5);
    4: two blocks with surplus proper structures (R = 5). See the appendix of the article for the
    generating procedure.
    """
    if case not in _SIM:
        raise ValueError("case must be 1, 2, 3 or 4")
    blocks = [_read(f"sim{case}_block{k}.csv", header=False)[2] for k in range(1, _SIM[case] + 1)]
    return dict(blocks=blocks, R=_SIM_R[case], block_names=[f"block {k}" for k in range(1, _SIM[case] + 1)])
