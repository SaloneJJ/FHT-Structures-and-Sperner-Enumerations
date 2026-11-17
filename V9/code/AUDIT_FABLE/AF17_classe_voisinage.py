#!/usr/bin/env python3
"""AF17 -- Elargissement exhaustif (au sens defini) de la classe K26.  Criteres et definition de « exhaustif » : 05_audit/UNICITE_ET_CLASSE_14-11.md,
section 3, ecrits le 14/11/2026 a 14 h 50 AVANT calcul.

(1) bfs : tous les membres de la classe a distance <= R (mouvements d'orbite : suppression, reduction, agrandissement elementaire, comme AF15/A2)
    des representants connus (K^(4), K26°, K26°°, K26*, V129-011/-017/-018, + fichiers donnes avec --depart).  Appartenance : recouvrement,
    Betti Z/2, cone de l'inclusion acyclique, |Aut| = 13 824 (A2, Prop. « Orbital moves keeping the class »).  Exhaustif a distance R.
    Etat sauve apres chaque sommet developpe (reprise automatique) : AF17_bfs_etat.json ; membres : AF17_membres.jsonl.
(2) deco : decorations a un niveau de chaque representant (ou de chaque membre trouve, option --membres) :
    pendant (nouveau sommet n, facette {n, u}) et jumeau (u' ajoute a toutes les facettes contenant u), sur un sommet u representant
    de chaque orbite de sommets, ou sur toute l'orbite ; d'un seul cote, ou des deux cotes (image par sigma).  Test complet AF9 (Aut recalcule).
    Sortie : AF17_deco.jsonl.
Usage : python3 AF17_classe_voisinage.py bfs [R=1] [--depart f.json ...]
        python3 AF17_classe_voisinage.py deco [--membres]
        python3 AF17_classe_voisinage.py bilan
"""
import os, sys, json, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import AF15_paysage_puits as P15
from AF15_paysage_puits import load, moves, apply, in_class, names, tomask, V, n, betti_z2, Key, fmt
CAT = os.path.join(HERE, '..', '..', 'catalogues')
ETAT = os.path.join(HERE, 'AF17_bfs_etat.json'); MEMB = os.path.join(HERE, 'AF17_membres.jsonl'); DECO = os.path.join(HERE, 'AF17_deco.jsonl')

def starts(extra):
    S = {s: load(s) for s in ('K4', 'V36', 'V37', 'V32')}
    for l in open(os.path.join(CAT, 'V129_candidats_facettes.jsonl')):
        d = json.loads(l)
        if d['id'] in ('V129-011', 'V129-017', 'V129-018'): S[d['id']] = [tomask(f) for f in d['facets']]
    for f in extra:
        d = json.load(open(f)); d = d['facettes'] if isinstance(d, dict) and 'facettes' in d else (d['facets'] if isinstance(d, dict) else d)
        S[os.path.basename(f)] = [tomask(x) for x in d]
    return S

def key(F): return ','.join(map(str, sorted(F)))

def bfs(R, extra):
    if os.path.exists(ETAT):
        st = json.load(open(ETAT)); print(f'reprise : {len(st["vus"])} membres vus, file {len(st["file"])}', flush=True)
    else:
        S = starts(extra); st = {'vus': {}, 'file': []}
        for nm, F in S.items():
            k = key(F)
            if k not in st['vus']:
                st['vus'][k] = {'nom': nm, 'dist': 0, 'parent': None}; st['file'].append(k)
                with open(MEMB, 'a') as fo: fo.write(json.dumps({'nom': nm, 'dist': 0, 'parent': None, 'facets': [names(f) for f in F]}) + '\n')
    b0 = {3: 12, 4: 36, 5: 36, 6: 12}
    t0 = time.time(); cpt = collections.Counter(d['dist'] for d in st['vus'].values())
    while st['file']:
        k = st['file'][0]; info = st['vus'][k]
        if info['dist'] >= R: st['file'].pop(0); continue
        F = [int(x) for x in k.split(',')]; t = time.time(); new = 0; tally = collections.Counter()
        for kind, o in moves(F):
            nE, small = apply(F, kind, o)
            if nE is None: tally['non elementaire'] += 1; continue
            kk = key(nE)
            if kk in st['vus']: tally['deja vu'] += 1; continue
            ver = in_class(F, nE, small, b0); tally[ver] += 1
            if ver != 'classe': continue
            nm = f'M{len(st["vus"])}'; st['vus'][kk] = {'nom': nm, 'dist': info['dist'] + 1, 'parent': info['nom']}; st['file'].append(kk); new += 1
            with open(MEMB, 'a') as fo:
                fo.write(json.dumps({'nom': nm, 'dist': info['dist'] + 1, 'parent': info['nom'], 'mouvement': kind,
                                     'kappa': len(nE), 'volume': P15.volume(nE), 'facets': [names(f) for f in nE]}) + '\n')
        st['file'].pop(0)
        json.dump(st, open(ETAT + '.tmp', 'w')); os.replace(ETAT + '.tmp', ETAT)
        cpt = collections.Counter(d['dist'] for d in st['vus'].values())
        print(f"{info['nom']} (d={info['dist']}) : {new} nouveaux membres ; {dict(tally)} ; total par distance {dict(sorted(cpt.items()))} ; "
              f"file {len(st['file'])} ({time.time()-t:.0f} s, cumul {(time.time()-t0)/3600:.2f} h)", flush=True)
    print(f'=== BFS termine (R = {R}) : membres par distance {dict(sorted(cpt.items()))}', flush=True)

def decorations(F):
    """(nom, facettes en noms) pour chaque decoration a un niveau"""
    Fn = [names(f) for f in F]
    vorb, seen = [], set()
    G = P15.GRP
    for i in range(n):
        if i in seen: continue
        o = sorted({p[i] for p in G}); seen |= set(o); vorb.append([V[j] for j in o])
    def sw(v): return v.replace('L_', '\x00').replace('R_', 'L_').replace('\x00', 'R_')
    out = []
    for orb in vorb:
        for portee in (('un', 'orbite') if len(orb) > 1 else ('un',)):
            targets = [orb[0]] if portee == 'un' else orb
            for cote in ('un_cote', 'deux_cotes'):
                if cote == 'deux_cotes':
                    if all(sw(u) == u for u in targets): continue
                    T = targets + [sw(u) for u in targets]
                else: T = targets
                for typ in ('pendant', 'jumeau'):
                    nF = [list(f) for f in Fn]
                    for u in T:
                        if typ == 'pendant': nF.append(['P_' + u, u])
                        else: nF = [f + ['T_' + u] if u in f else f for f in nF]
                    out.append((f'{typ}:{portee}:{cote}:{orb[0]}', nF))
    return out

def deco(use_members):
    from AF9_porteurs_triplet import test, verdict
    from AF3_K26_homologie_automorphismes import K26_definition
    base = starts([])
    if use_members and os.path.exists(MEMB):
        for l in open(MEMB):
            d = json.loads(l); base.setdefault(d['nom'], [tomask(f) for f in d['facets']])
    done = set()
    if os.path.exists(DECO):
        for l in open(DECO): done.add(json.loads(l)['id'])
    t0 = time.time(); o, ref = test('K26', K26_definition(), want_chars=True)
    with open(DECO, 'a') as fo:
        for nm, F in base.items():
            for dn, nF in decorations(F):
                ident = f'{nm}|{dn}'
                if ident in done: continue
                t = time.time(); o, _ = test(ident, nF, ref_chars=ref); ver = verdict(o)
                o = {k: (list(v) if isinstance(v, (set, tuple)) else v) for k, v in o.items()}; o.update(id=ident, verdict=ver)
                fo.write(json.dumps(o, default=str) + '\n'); fo.flush()
                print(f'{ident:45s} v={o["v"]} |Aut|={o["aut"]} -> {ver} ({time.time()-t:.0f} s)', flush=True)
    print(f'fin decorations ({(time.time()-t0)/3600:.2f} h)')

def bilan():
    if os.path.exists(MEMB):
        M = [json.loads(l) for l in open(MEMB)]
        c = collections.Counter(d['dist'] for d in M); print('membres (BFS) par distance :', dict(sorted(c.items())))
        sm = sorted((d for d in M if 'volume' in d), key=lambda d: d['volume'])[:10]
        print('plus petits volumes :', [(d['nom'], d['volume'], d['kappa']) for d in sm])
    if os.path.exists(DECO):
        D = [json.loads(l) for l in open(DECO)]
        c = collections.Counter(d['verdict'] for d in D); print('decorations :', dict(c))
        print('decorations porteuses :', [d['id'] for d in D if d['verdict'].startswith('PORTEUR')])

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: print(__doc__); sys.exit()
    extra = a[a.index('--depart') + 1:] if '--depart' in a else []
    if a[0] == 'bfs': bfs(int(a[1]) if len(a) > 1 and a[1].isdigit() else 1, extra)
    elif a[0] == 'deco': deco('--membres' in a)
    elif a[0] == 'bilan': bilan()
