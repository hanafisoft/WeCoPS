"""Reproduce the 28 verdicts of the article (four weightings, seven cases).

Run:  python reproduce_article.py
"""
import wecops as wc

CASES = [
    ("Simulated 1", wc.load_simulated(1)["blocks"], 5),
    ("Simulated 2", wc.load_simulated(2)["blocks"], 10),
    ("Simulated 3", wc.load_simulated(3)["blocks"], 5),
    ("Simulated 4", wc.load_simulated(4)["blocks"], 5),
    ("Wines", wc.load_wines()["blocks"], 3),
    ("Butters, R = 4", wc.load_butters()["blocks"], 4),
    ("Butters, R = 3", wc.load_butters()["blocks"], 3),
]
WEIGHTINGS = ["rv", "mfa", "statis", "nvar"]

n_eq = n_fa = 0
print(f"{'case':16s}" + "".join(f"{w:>14s}" for w in WEIGHTINGS))
for name, blocks, R in CASES:
    cells = []
    for w in WEIGHTINGS:
        r = wc.compare_weighting(blocks, R, w)
        fa = r.equivalent and not r.pc_same
        n_eq += r.equivalent
        n_fa += fa
        cells.append(("=" if r.equivalent else "#") + f" {r.kappa.min():.2f}" + (" c" if fa else "  "))
    print(f"{name:16s}" + "".join(f"{c:>14s}" for c in cells))
print(f"\nequivalent: {n_eq} / 28; false alarms of the principal-component reading: {n_fa} / {n_eq}")
print("(= equivalent, # interpretation changed, value: min presence, c: false alarm)")
