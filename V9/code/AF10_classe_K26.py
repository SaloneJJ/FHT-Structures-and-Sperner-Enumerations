#!/usr/bin/env python3
"""AUDIT -- Classe K26 elargie (decision de JJ du 14/11/2026) : caracteristiques propres de chaque membre.
Membres : representants (memes 26 sommets, |Aut| = 13 824, meme triplet : AF9) et decores (pendant, jumeau : homologie inchangee,
groupe de generation eventuellement reduit).  Pour chacun : sommets, aretes (facettes), simplexes, |Aut|, ordre du groupe induit sur
G0 = {g0..g3} (groupe de generation), orbites de facettes, taille maximale d'une arete, rayon et diametre du graphe primal,
nombre de voisins de L_z et de R_z (sommet « central » ou non), sommets hors des 26 sommets standard (decoration).
Python 3 + numpy ; reutilise AF9 (automorphismes) et AF3 (Definition).  N'ecrit aucun fichier.
Usage : python3 AF10_classe_K26.py RACINE_FHT26
"""
import json, os, sys
from collections import deque
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import AF9_porteurs_triplet as T
from AF3_K26_homologie_automorphismes import K26_definition

STD = set(T.Cx('K26', K26_definition()).V)

def primal(C):
    nb = {v: set() for v in range(len(C.V))}
    for f in C.F:
        for a in f:
            nb[a] |= set(f) - {a}
    return nb

def ecc(nb, s):
    d = {s: 0}; q = deque([s])
    while q:
        u = q.popleft()
        for w in nb[u]:
            if w not in d: d[w] = d[u] + 1; q.append(w)
    return max(d.values())

def facet_orbits(C, gens):
    F = list(C.F); ix = {f: i for i, f in enumerate(F)}; par = list(range(len(F)))
    def find(i):
        while par[i] != i: par[i] = par[par[i]]; i = par[i]
        return i
    for g in gens:
        for i, f in enumerate(F):
            j = ix[frozenset(g[v] for v in f)]
            a, b = find(i), find(j)
            if a != b: par[a] = b
    return len({find(i) for i in range(len(F))})

def describe(name, facets, role):
    C = T.Cx(name, facets); gens, order = T.aut_sgs(C); nv = len(C.V)
    G = T.enumerate_group(gens, nv)
    g0 = [C.ix[v] for v in ('g0', 'g1', 'g2', 'g3') if v in C.ix]
    gen_img = len({tuple(g[i] for i in g0) for g in G}) if G is not None else None
    nb = primal(C); e = [ecc(nb, v) for v in range(nv)]
    deg = lambda lab: len(nb[C.ix[lab]]) if lab in C.ix else None
    extra = sorted(set(C.V) - STD)
    return dict(name=name, role=role, v=nv, aretes=len(C.F), simplexes=C.nsimp, aut=order, gen=gen_img,
                orb=facet_orbits(C, gens), fmax=max(len(f) for f in C.F), rayon=min(e), diam=max(e),
                Lz=deg('L_z'), Rz=deg('R_z'), extra=extra)

if __name__ == '__main__':
    root = sys.argv[1]; cat = os.path.join(root, '03_code', 'catalogues'); obj = os.path.join(root, '04_objets')
    def load(fn):
        d = json.load(open(os.path.join(obj, fn))); return d if isinstance(d, list) else d.get('facets', d.get('facettes'))
    et = json.load(open(os.path.join(cat, 'V145_etoiles_facettes.json')))['etoiles']
    v129 = {json.loads(l)['id']: json.loads(l)['facets'] for l in open(os.path.join(cat, 'V129_candidats_facettes.jsonl'))}
    members = [('K^(4)', K26_definition(), 'representant (Definition de A2)'),
               ('K26° (V36)', load('V36_FHT_K26vol.json'), 'representant (descente en volume)'),
               ('K26°° (V37)', load('V37_FHT_K26volmin.json'), 'representant (recuit ; plus petit volume trouve)'),
               ('K26* (V32)', load('V32_FHT_K26etoile.json'), 'representant (descente en aretes)'),
               ('V129-011', v129['V129-011'], 'representant (V129 : ajout/echange depuis K^(4))'),
               ('V129-017', v129['V129-017'], 'representant (V129 : agrandissement depuis K26°°)'),
               ('V129-018', v129['V129-018'], 'representant (V129 : agrandissement depuis V129-017)'),
               ('K_27_546', et['K_27_546']['facets'], 'decore (pendant)'),
               ('K_36_44', et['K_36_44']['facets'], 'decore (pendant + jumeau ; sigma par prefixe)')]
    print('| membre | role | sommets | aretes | simplexes | |Aut| | groupe sur G0 | orbites d\'aretes | arete max | rayon / diam. | voisins L_z / R_z | sommets ajoutes |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|')
    for nm, fac, role in members:
        d = describe(nm, fac, role)
        print(f"| {d['name']} | {d['role']} | {d['v']} | {d['aretes']} | {d['simplexes']} | {d['aut']} | {d['gen']} | {d['orb']} | {d['fmax']} | "
              f"{d['rayon']} / {d['diam']} | {d['Lz']} / {d['Rz']} | {', '.join(d['extra']) or '—'} |", flush=True)
