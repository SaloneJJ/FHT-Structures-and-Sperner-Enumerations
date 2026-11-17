# Trois vérifications avant la soumission EDSU-2027 — critères écrits AVANT tout calcul
16 nov. 2026, ~18 h 45 UTC. Aucun calcul n'a été fait avant la rédaction de ce fichier.

## Données (relevées en ligne le 16/11/2026 ; les tables Δχ² de nu-fit.org ne sont pas accessibles depuis cette session)
- NuFIT 6.1 (données jusqu'à nov. 2025, JUNO inclus pour θ12 et Δm²21), variante IC24 + SK-atm, telle que reprise dans arXiv:2608.18812 et arXiv:2609.01196 :
  - NO : sin²θ23 = 0.470 (+0.017 −0.013), 3σ [0.435, 0.585] ; δ = 212° (+26 −41), 3σ [124°, 364°].
  - IO : sin²θ23 = 0.550 (+0.012 −0.015), 3σ [0.440, 0.584] ; δ = 274° (+22 −25), 3σ [201°, 335°].
  - sin²θ12 = 0.3088 (+0.0067 −0.0066) ; sin²θ13 = 0.02248 (NO), 0.02262 (IO).
  - Δm²31 = 2.511 (+0.021 −0.020) ×10⁻³ eV² (NO, avec SK) ; Δm²32 = −2.484 ± 0.020 (IO).
  - Δm²21 : 7.530 (+0.096 −0.097) ×10⁻⁵ eV² (JUNO + NuFIT 6.1, arXiv:2601.09791).
- JUNO, 59,1 jours : sin²θ12 = 0.3092 ± 0.0087.
- Capozzi et al. (PRD 2026), combinaison globale : sin²θ12 = 0.3085 ± 0.0073.
- Ordre des masses : avec JUNO, l'ordre inversé est défavorisé de Δχ² = 4.6 (sans SK) / 9.4 (avec SK-atm + IC24) [2601.09791] ;
  DESI DR2 + CMB (ΛCDM) : Σmν < 0.0642 eV (95 %), contre un minimum de 0.099 eV en IO.

## Méthode déclarée (approximation)
Faute des tables, Δχ² à une dimension est **estimé** par interpolation linéaire de √Δχ² entre le meilleur ajustement (0), le bord à 1σ (1) et le bord à 3σ (3), de chaque côté.
Au-delà du bord à 3σ : « Δχ² > 9 », sans extrapolation.
Limites :
- les profils réels ne sont pas paraboliques ;
- le profil de θ23 en NO est bimodal (minimum local vers 0.56), donc la valeur estimée à 0.5 est incertaine d'environ un facteur 2 ;
- les corrélations θ23–δ ne sont pas prises en compte, et l'on ne somme pas les Δχ² à une dimension.
Le résultat est à confirmer avec les tables Δχ² 1D et 2D de NuFIT 6.1 (script à fournir quand JJ les aura téléchargées).

## V1 — θ23 et δ_CP
- Points à tester :
  - TM1 exact : sin²θ23 = 0.5 et δ = 90° ou 270° (signe non prédit : on retient le meilleur des deux) ;
  - Georgi–Jarlskog (GJ) : sin²θ23 ∈ [0.493, 0.510] (min du Δχ² sur l'intervalle) et |δ| sur la courbe de la Prop. corr au sin²θ12 mesuré.
- Verdict par paramètre et par ordre : Δχ² ≤ 4 « compatible » ; 4 < Δχ² ≤ 9 « en tension » ; Δχ² > 9 « défavorisé ».
- Conséquence pour A2 :
  - si δ = 90° est défavorisé et δ = 270° compatible dans un ordre, A2 doit dire que les données choisissent δ = −90°, ce qui fixe l'orientation (quel vecteur propre est le muon) ;
  - si les deux signes sont défavorisés dans les deux ordres, les deux scénarios sont en difficulté, et A2 doit le dire.

## V2 — ordre des masses et m_ββ en ordre normal
- Recalculer avec les entrées NuFIT 6.1 :
  - les bandes de la Prop. mbb (ii) en NO, pour m1 → 0, sur la boîte des variations indépendantes à 3σ de Δm²21, Δm²31 et s13² ;
  - la valeur (i), avec |Ue1|² = c12² c13² et sin²θ12 dans son intervalle à 3σ ;
  - Σmν minimal en NO et en IO.
- Critère : si un bord de bande bouge de plus de 0.2 meV, mettre A2 à jour ; sinon, le dire.
- Conséquence : l'ordre inversé est défavorisé par l'oscillation (≈ 2–3σ) et par la cosmologie (ΛCDM). La phrase d'A2 sur l'ordre est actualisée avec ces chiffres, et les prédictions NO deviennent l'énoncé principal pour m_ββ.

## V3 — sin²θ12
- Prédiction TM1 : sin²θ12 = 1 − 2/(3 c13²), avec s13² de NuFIT 6.1 (NO) et sa barre d'erreur.
- Écart en σ (erreur de mesure et erreur de la prédiction en quadrature) face à JUNO, à NuFIT 6.1 et à Capozzi et al.
  Verdict : < 1σ compatible ; 1–2σ « au bord » ; 2–3σ tension ; > 3σ exclu.
- GJ : |δ| lu sur la courbe pour sin²θ12 à ±1σ de NuFIT 6.1 (interpolation linéaire de la table de la Prop. corr).

## Addendum (16 nov., ~19 h 15 UTC), écrit après réception des tables NuFIT 6.1 et avant leur lecture numérique
Seuls les en-têtes des fichiers et la notice (v61.release-notes.pdf) ont été lus.
- Les tables officielles remplacent l'interpolation de V1 : c'est un changement de méthode, déclaré ici. Les seuils restent inchangés.
- Identification des fichiers (téléchargés sans nom) : par l'en-tête « ### Normal/Inverted Ordering ». Pour distinguer TByes (avec SK-atm + IC24) de TBoff (IC23, sans), on utilise deux signatures connues :
  - le minimum de Δχ² du fichier IO, plus grand avec TByes ;
  - la position du minimum de sin²θ23 en NO, attendue vers 0.47 avec TByes et vers 0.56 sans.
  Si ces deux signatures ne concordent pas, on s'arrête et on le signale.
- Δχ² = valeur de la table moins le minimum global sur les deux ordres (convention NuFIT). On rapporte aussi la valeur relative au minimum de l'ordre considéré.
- Lectures :
  - 1D T23 en 0.5 et min sur [0.493, 0.510] ;
  - 1D DCP en ±90° (la table utilise δ ∈ [−180°, 180°]) ;
  - 1D T12 en 0.318 ;
  - **2D T23/DCP** en (0.5, ±90°) et min sur le rectangle GJ [0.493, 0.510] × ±[89°, 98°] ;
  - **2D T12/DCP** le long de la courbe GJ (|δ| ↔ sin²θ12, table de la Prop. corr, signe libre) : min du Δχ².
  - Interpolation bilinéaire sur la grille, sans lissage.
- Verdicts :
  - 1D : ≤ 4 compatible, ≤ 9 en tension, > 9 défavorisé ;
  - 2D (2 degrés de liberté) : ≤ 6.18 compatible (2σ), ≤ 11.83 en tension (3σ), au-delà défavorisé.
- On ne somme pas de projections différentes.

## Addendum 2 (~19 h 25 UTC) : identification des fichiers NO
La signature « argmin sin²θ23 » ne discrimine pas : les deux fichiers NO donnent 0.470, et le script a attribué le même fichier aux deux variantes, ce qui est une erreur de code.
Les résultats déjà imprimés pour NO sont donc à refaire. La signature IO, elle, discrimine (minima 1.49 et 5.91).
Nouvelle signature, déclarée ici avant lecture : le meilleur ajustement de Δm²3ℓ publié pour NuFIT 6.1 dans arXiv:2601.09791 (éq. 13 et 13a) :
- NO : 2.521 sans SK-atm, 2.511 avec ;
- IO : −2.510 sans, −2.484 avec.
Un fichier est attribué à la variante dont la valeur est la plus proche, à condition que l'écart entre les deux candidats dépasse 0.005 ; sinon, on s'arrête.
