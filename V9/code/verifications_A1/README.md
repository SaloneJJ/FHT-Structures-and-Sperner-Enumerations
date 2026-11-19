# Verifications for Article I (Fractal Hyper-Trees I. Combinatorial foundations), revised version

Independent checks added during the revision of Article I (November 2026). Python 3 with numpy and scipy.
Each script prints its results; the `.log` files are the outputs of the runs cited in the article.

| script | what it checks (Article I) | run |
|---|---|---|
| `test_axioms_literal.py` | literal test of the axioms (A1)–(A5) over all edge families and all level functions with values in {1,…,v−1}, v ≤ 4: 1, 4, 23 isomorphism classes (Section "Enumeration", validation) | `python3 test_axioms_literal.py` (a few seconds) |
| `twisted_decomp.py` | numerical Wedderburn decomposition of real twisted group algebras R^η[G] (centre, blocks of a generic central element, type R/C/H by the signature of the trace form); checks on R^sgn[S4] = A_F and R^sgn[S5]; decompositions for S4×S3 and S4×C2 and all their sign characters (census of S4, sectors) | `python3 twisted_decomp.py` |
| `valpha54_check.py` | recomputes, for the representative groups of the 54 classes (v,α), v ≤ 6 (generators read from `valpha54_table_rows.tex`, extracted from the table of the article), the type (r,c,h) of R[G] and the decomposition of R^η[G] for every sign character; lists the sectors A_F | `python3 valpha54_check.py` (needs `twisted_decomp.py`) |
| `gl23_check.py` | R^det[GL(2,3)] = A_F ⊕ H ⊕ H ⊕ M2(H); det = sign of the permutation of the 9 points of F_3^2; counts of g² = ±I used in the proof of the spinorial sector | `python3 gl23_check.py` (needs `twisted_decomp.py`) |
| `norms_check.py` | character norms: Λ² of the standard representation of S5 and spin 3/2 restricted to the binary octahedral group 2O are irreducible; degrees 2, 2, 4 of the faithful representations | `python3 norms_check.py` |

These computations are in floating point; they are checks, not exact certifications.
Written with the assistance of Claude (Anthropic), checked and run by the author. Code: MIT; logs: CC BY 4.0.
