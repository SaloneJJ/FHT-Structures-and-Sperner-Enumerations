import sys, os
D = sys.argv[1]
for f in sorted(os.listdir(D)):
    rows, on, order = [], False, None
    for line in open(os.path.join(D, f)):
        if line.startswith('###'): order = line.strip('# \n')
        if line.startswith('# DMA projection'): on = True; continue
        if on and line.startswith('#'): break
        if on and line.strip(): rows.append([float(x) for x in line.split()])
    r = min(rows, key=lambda t: t[1]); print(f, order, 'argmin Δm²3ℓ =', r[0], 'Δχ² =', r[1])
