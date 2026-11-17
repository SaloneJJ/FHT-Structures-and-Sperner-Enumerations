# Résultats des trois vérifications (critères : CRITERES_VERIF_EDSU.md ; script verif_edsu.py ; journal verif_edsu.log)
16 nov. 2026, ~19 h UTC. Changement déclaré après le premier calcul : ajout de « V2 bis » (m_ββ en NO avec m1 borné par DESI), qui ne modifie aucun verdict.

## V1 — θ23 et δ (Δχ² 1D estimés depuis les intervalles 1σ/3σ publiés de NuFIT 6.1, IC24 + SK-atm)
| | sin²θ23 = 0.5 | GJ [0.493, 0.510] | δ = +90° | δ = −90° (270°) | GJ, |δ| ∈ [89°, 99°], signe libre |
|---|---|---|---|---|---|
| NO | ≈ 1.6 compatible | ≈ 1.3 | > 9 défavorisé | ≈ 2.3 compatible | ≈ 1.9 (δ = 261°) |
| IO | ≈ 3.0 compatible | ≈ 2.3 | > 9 défavorisé | ≈ 0.0 compatible | ≈ 0.0 (δ = 271°) |
Ordre : IO défavorisé, Δχ² = 4.6 / 9.4 (JUNO + NuFIT 6.1).
**Verdict** : les deux scénarios sont compatibles, à condition que δ = −90°.
Le signe, non prédit, est donc fixé par les données : il dépend de l'orientation (lequel de v_ω, v_ω̄ est le muon).
Réserve : approximation, et profil de θ23 bimodal en NO. À confirmer avec les tables Δχ² de nu-fit.org, qui ne sont pas accessibles d'ici (robots, liste blanche réseau).

## V2 — ordre et m_ββ
- NO, m1 → 0, TM1 exact (ii) : 1.40–1.74 ou 3.69–3.96 meV (A2 V9 : 1.4–1.8 ou 3.6–4.0).
  Écarts < 0.2 meV, donc la Prop. reste inchangée ; les tableaux citent maintenant les valeurs NuFIT 6.1.
- Réalité seule (i), m1 → 0 : 1.17–4.07 meV (A2 : 1.0–4.2), là aussi écarts < 0.2 meV.
- Σmν minimal : NO 58.8 meV, IO 98.9 meV. DESI DR2 + CMB (ΛCDM) : < 64.2 meV, d'où m1 < 4.25 meV en NO.
- V2 bis : pour m1 ∈ [0, 4.25] meV, les quatre courbes TM1 restent dans [0, 6.97] meV, et l'une d'elles atteint 0.
**Verdict** : l'énoncé principal sur m_ββ devient celui de l'ordre normal (< 7 meV, hors de portée de la prochaine génération) ; les bandes de l'ordre inversé passent au conditionnel.

## V3 — sin²θ12
- TM1 : 0.3180 ± 0.0004. Écarts : +1.01σ (JUNO), +1.37σ (NuFIT 6.1), +1.30σ (Capozzi et al.), +1.47σ (JUNO + NuFIT 6.1).
  Verdict : « au bord » (1–2σ), en accord avec Zhang, PRD 113 (2026) 055035, pour qui TM1 est au bord de la région à 1σ.
- GJ : à ±1σ de NuFIT 6.1, |δ| ≈ 89°–98° ; au centre, 92°–95°.

## Modifications d'A2 V9 (copie avant modification : A2_FiniteTriple/_archive_V9_avant_verif_EDSU/)
- Résumé (5) : δ = −90° choisi par les données ; θ23 = 45° à ≈ 2σ ; m_ββ en NO en premier (< 7 meV ; 1.4–1.7 ou 3.7–4.0 meV), bandes IO au conditionnel.
- Tables des scénarios et des prédictions : colonnes « données » (NuFIT 6.1, JUNO, écarts en σ, signe de δ, m_ββ en NO).
- Paragraphe « What is firm » : ordre (Δχ² 4.6 / 9.4, DESI 0.064 eV) ; statut face aux données, avec l'approximation déclarée.
- Paragraphe « Tests » : cosmologie (DESI DR2, relâchée si l'énergie noire évolue) ; le NO hors de portée de la prochaine génération 0νββ.
- Bibliographie : NuFit61, DESI25nu, Zhang26.
- Compilation : 71 pages, 0 référence indéfinie, 11 overfull (inchangé).

---
# Mise à jour (~19 h 35 UTC) : tables officielles NuFIT 6.1 (script verif_nufit61.py ; journal verif_nufit61.log)
Fichiers (nufit61/, noms de téléchargement conservés), identifiés par Δm²3ℓ (addendum 2) :
- TBoff = 144289 (NO), 196693 (IO) ;
- TByes = 506694 (NO), 460581 (IO) ;
- 838518 = notice.
Incident déclaré : la première identification (argmin θ23) ne discriminait pas les fichiers NO, et le premier passage a lu deux fois le même fichier NO. Ce passage est annulé et refait.
Ordre des masses dans NuFIT 6.1 seul : Δχ²(IO) = 1.49 (TBoff), 5.91 (TByes).

## Valeurs, TByes (variante principale). Δχ² relatifs au minimum global ; [ ] = relatif au minimum IO
| | NO | IO |
|---|---|---|
| 1D sin²θ23 = 0.5 | 2.12 | 8.88 [2.97] |
| 1D δ = +90° | 14.79 | 59.5 [53.6] |
| 1D δ = −90° | 2.72 | 5.93 [0.03] |
| 1D sin²θ12 = 0.318 | 1.84 | [1.85] |
| 2D (θ23, δ) = (0.5, −90°), TM1 exact | **7.12** (≈ 2.2σ, 2 ddl) — en tension | 8.88 [2.97] |
| 2D GJ (θ23 ∈ [0.493, 0.510], \|δ\| ∈ [89°, 98°]) | **5.56** (≈ 1.9σ) — compatible | 8.59 [2.68] |
| 2D (θ12, δ) = (0.318, −90°) | 4.55 | 7.78 [1.87] |
| 2D (θ12, δ), courbe GJ | 2.51 en (0.308, −95°) | 5.98 [0.07] |
TBoff, NO : (0.5, −90°) = 8.31 ; GJ = 6.75 ; δ = +90° : 11.2 ; δ = −90° : 3.7.

## Ce qui change par rapport à l'estimation
- L'interpolation sous-estimait la tension conjointe. Séparément, θ23 = 45° (2.1) et δ = −90° (2.7) sont compatibles ; ensemble, en NO, la corrélation porte le point TM1 exact à Δχ² = 7.1 (2.2σ).
- Verdict : les deux scénarios survivent à environ 2σ en NO (GJ à 1.9σ, TM1 exact à 2.2σ). δ = +90° est exclu dans les deux ordres. En IO, à l'intérieur de cet ordre, les deux sont très compatibles, mais l'IO est lui-même défavorisé.
- A2 V9 corrigé en conséquence (résumé, deux tables, paragraphe « Status against the data ») ; 71 pages, 0 référence indéfinie.

## Contrôle d'identification (~19 h 40 UTC)
JJ a redéposé les fichiers sous leurs noms officiels (« v61.release-… (1).txt.xz »).
Contenu décompressé comparé par sha256 aux fichiers utilisés :
- TBoff-NO = 144289, TBoff-IO = 196693, TByes-NO = 506694, TByes-IO = 460581 : tous identiques.
L'identification par Δm²3ℓ est donc confirmée, et les résultats ci-dessus sont inchangés.
