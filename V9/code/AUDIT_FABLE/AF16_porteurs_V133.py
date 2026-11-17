#!/usr/bin/env python3
"""AF16 -- Test complet des porteurs du triplet sur la table d'exploration V133 (FHT_22).
Criteres : ceux de AF9 (05_audit/TEST_PORTEURS_TRIPLET.md, ecrits le 14/11/2026 a 7 h 05) ; liste des candidats fixee le 14/11 a 14 h 45 AVANT calcul
(05_audit/UNICITE_ET_CLASSE_14-11.md, section 2).
Selection dans la table V52_paysage_v5_f30.json (31 394 enregistrements au 11/11, 21 h 43) :
   s4cube : tests.gamma = ok, tests.J_nec = ok et |Aut| = 13 824          (123 ; les « 103 » de V133 au 11/11 matin)
   jnec   : tests.gamma = ok et tests.J_nec = ok, tout |Aut|                (1 598)
   reste  : jnec sans les 123 de s4cube                                     (1 475) -- a lancer plutot que jnec
Pour chacun : AF9.test (sigma, Aut, homologie sur F_65521, J' exact C1, caracteres C2 compares a K^(4)) -> verdict.
Reprise automatique : les objets deja dans le fichier de sortie sont sautes.  Decoupage : --part i/n (i = 0..n-1) pour plusieurs terminaux.
Usage : python3 AF16_porteurs_V133.py s4cube|reste|jnec [--part i/n] [--fht22 CHEMIN/FHT_22]
Sorties : AF16_<mode>[_part].jsonl (une ligne par objet) ; resume : python3 AF16_porteurs_V133.py bilan s4cube|jnec
"""
import os, sys, json, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
args = sys.argv[1:]
def opt(name, default=None):
    if name in args:
        i = args.index(name); v = args[i + 1]; del args[i:i + 2]; return v
    return default
FHT22 = opt('--fht22', os.path.normpath(os.path.join(HERE, '..', '..', '..', '..', 'FHT_22')))
PART = opt('--part')
CODE = os.path.join(FHT22, '05_code', 'code_v30')
TAB, JL, IDX = (os.path.join(CODE, f) for f in ('V52_paysage_v5_f30.json', 'V133_FHT_facettes.jsonl', 'V133_FHT_facettes.jsonl.idx'))

def candidates(mode):
    R = json.load(open(TAB))['recs']
    out = []
    for r in R:
        t = r.get('tests') or {}
        if not str(t.get('gamma', '')).startswith('ok') or not str(t.get('J_nec', '')).startswith('ok'): continue
        if mode == 's4cube' and r['aut'] != 13824: continue
        if mode == 'reste' and r['aut'] == 13824: continue      # les 1 475 restants, sans refaire les 123 de s4cube
        out.append((r['id'], r['v'], r['aut'], r['vol']))
    return sorted(out, key=lambda x: (x[3], x[0]))          # du plus petit volume au plus grand

def facets_of(ident, off, fh):
    fh.seek(off); d = json.loads(fh.readline())
    assert d['id'] == ident, (d['id'], ident)
    return d['facets']

def outname(mode, part):
    return os.path.join(HERE, f'AF16_{mode}' + (f'_part{part.replace("/", "de")}' if part else '') + '.jsonl')

def run(mode):
    from AF9_porteurs_triplet import test, verdict
    from AF3_K26_homologie_automorphismes import K26_definition
    C = candidates(mode)
    if PART:
        i, n = map(int, PART.split('/')); C = C[i::n]
    out = outname(mode, PART)
    done = set()
    if os.path.exists(out):
        for l in open(out):
            try: done.add(json.loads(l)['id'])
            except Exception: pass
    print(f'{len(C)} candidats ({mode}{", part " + PART if PART else ""}) ; deja faits : {len(done & {c[0] for c in C})}', flush=True)
    t0 = time.time(); o, ref = test('K26', K26_definition(), want_chars=True)
    print(f'reference K^(4) : |Aut| {o["aut"]}, Betti {o["betti"]}, C1 {o["C1"]} ({time.time()-t0:.0f} s)', flush=True)
    off = json.load(open(IDX))['off']
    with open(JL) as fh, open(out, 'a') as fo:
        for k, (ident, v, aut, vol) in enumerate(C):
            if ident in done: continue
            if ident not in off: print(f'{ident} : facettes absentes de l index', flush=True); continue
            t = time.time()
            try:
                o, _ = test(ident, facets_of(ident, off[ident], fh), ref_chars=ref); ver = verdict(o)
            except Exception as e:
                o, ver = {'name': ident}, f'ERREUR {type(e).__name__}: {e}'
            o = {k2: (v2 if not isinstance(v2, (set, tuple)) else list(v2)) for k2, v2 in o.items()}
            o.update(id=ident, verdict=ver, vol_table=vol)
            fo.write(json.dumps(o, default=str) + '\n'); fo.flush()
            print(f'[{k+1}/{len(C)}] {ident:10s} v={v} |Aut|={aut} vol={vol} -> {ver} ({time.time()-t:.0f} s)', flush=True)
    print(f'fin ({(time.time()-t0)/3600:.2f} h)', flush=True)

def bilan(mode):
    rows = {}
    for f in os.listdir(HERE):
        if f.startswith(f'AF16_{mode}') and f.endswith('.jsonl'):
            for l in open(os.path.join(HERE, f)):
                d = json.loads(l); rows[d['id']] = d
    for d in rows.values():      # requalification (14/11, 16 h 25) : C2 compare des ensembles d'etiquettes ; un groupe L/R reduit n'a pas toutes les etiquettes
        if d['verdict'].startswith('STRUCTURE REELLE') and d.get('C2_labels') and d['C2_labels'][0] < d['C2_labels'][1]:
            d['verdict'] = "J' EXACT MAIS GROUPE L/R REDUIT (A_F non lisible ; bimodule non compare)"
    c = collections.Counter(d['verdict'] for d in rows.values())
    print(f'{len(rows)} objets testes ({mode})')
    for k, n in c.most_common(): print(f'   {n:5d}  {k}')
    port = sorted((d for d in rows.values() if d['verdict'].startswith('PORTEUR')), key=lambda d: (d['v'], d['simplices']))
    print('porteurs (id, v, |Aut|, simplexes) :')
    for d in port: print(f'   {d["id"]:10s} {d["v"]:3d} {d["aut"]:7d} {d["simplices"]:7d}')
    autres = [d for d in rows.values() if d['verdict'].startswith('STRUCTURE REELLE')]
    if autres: print('structure reelle mais autre bimodule :', [d['id'] for d in autres])

if __name__ == '__main__':
    if not args: print(__doc__); sys.exit()
    if args[0] == 'bilan': bilan(args[1])
    elif args[0] == 'liste': C = candidates(args[1]); print(len(C)); print(collections.Counter(c[1] for c in C).most_common())
    else: run(args[0])
