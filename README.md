# Fractal Hyper-Trees (FHT): programs and logs

Programs, logs and verification scripts for the articles on Fractal Hyper-Trees by Jean-Jacques Salone
(CRREF, University of the Antilles; LIM, University of La Réunion), in particular
*Fractal Hyper-Trees I. Combinatorial foundations: definition, enumeration, invariants, morphisms and the algebra C ⊕ H ⊕ M3(C)*
(HAL hal-05748266, arXiv:2609.19895).

[![DOI](https://zenodo.org/badge/1350712518.svg)](https://doi.org/10.5281/zenodo.22163670)

Cite all versions: doi:10.5281/zenodo.22163670. Current release: v9.2.

## Where the code is
- `V9/code/`: the programs cited in the articles (version 9). `V9/README.md` describes the folder and `V9/INVENTAIRE.md` maps each name cited in an article to a file.
- `V9/code/FHT_Enumeration.py` and `V9/code/FHT_Enum_results.json`: the counts f(v) for v ≤ 9, with the counts by depth.
- `V9/code/verifications_A1/`: the checks added in the revision of Article I (moves and transition graphs, census of the FHTs with automorphism group containing S4 up to v = 9, separation signature, twisted group algebras). See its `README.md` for the command and the log of each check.

## Counts
f(v), v = 1..9: 1, 1, 4, 23, 231, 16947, 490083358, 1392195551804840146,
789204635842035050272964246705408655.
The values at v = 8, 9 use the Dedekind numbers D(7), D(8) and the terms a(8), a(9) of OEIS A003182.

## Running
Python 3.9 or later with numpy and scipy; pynauty for some automorphism computations (`V9/requirements.txt`).
Each script of `verifications_A1/` runs with `python3 <script>.py` and writes the log deposited beside it.

## Earlier definition (v4.0.0 and before)
`fht.py`, `requirements.txt` and `results/` at the root implement the definition of versions ≤ 4.0.0
(doi:10.5281/zenodo.22948863). They are kept for the record. Their counts (1, 1, 4, 21, 199, …) differ from the
current ones and are superseded by `V9/`. Versions ≤ 3.x over-counted I(v) for v ≥ 6.

## Licences
Code: MIT (`LICENSE`). Logs, data and documentation in `V9/`: CC BY 4.0 (`V9/LICENSE-docs.md`).
