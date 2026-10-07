#!/usr/bin/env python3
"""FINISHER-05x IEEG056 (Fusca, Siebenhuehner, Wang et al. 2023 Nat Commun; Dryad doi:10.5061/dryad.vdncjsxzn) -> derivative BIDS tree.

- sourcedata/dryad-vdncjsxzn/: all 4 original release files byte-for-byte (Dryad digest re-checked) + SHA256SUMS.
- sub-NN/ieeg/ (NN = row of the 68-patient SEEG arrays, 01..68):
    sub-NN_task-rest_desc-sync_connectivity.mat : that patient's PLV (complex, 24 freq x contacts x contacts), DFA,
        elec_dist, ref_mask, acceptcontHealthy, acceptcontEZ, parcel_assign, HealthyEZcriteria row, frequencies (values unchanged)
    sub-NN_task-rest_desc-dfa_contacts.tsv : DFA exponents, one row per contact, one column per frequency (full float repr)
    sub-NN_task-rest_contacts.tsv : per-contact parcel code/name and the non-EZ / EZ contact masks
- group/seeg/parcel_order.tsv, frequencies.tsv ; MEG_data.mat and Simulated_data.zip stay group-level in sourcedata/
  (no session->subject mapping is released for the 192 MEG sessions).
- Round-trip: every per-subject array re-read and compared exactly to the original cell; TSV values re-parsed == original.
Usage: f05x_build056.py <textdir>
"""
import json, os, shutil, sys
import numpy as np, scipy.io as sio
sys.path.insert(0, os.path.dirname(__file__))
from f05x_lib import *  # noqa

ID = 'IEEG056'; SLUG = 'dryad-vdncjsxzn'
SRC = f'{ROOT}/work/{ID}/sourcedata/{SLUG}'
TREE = f'{ROOT}/work/{ID}/bids'
REP = f'{ROOT}/work/{ID}/reports/f05x'
TEXT = sys.argv[1]
assert not os.path.exists(TREE + '/.nemar'), 'tree already uploaded; refusing to rebuild'
if os.path.exists(TREE):
    shutil.rmtree(TREE)
os.makedirs(TREE)
man = manifest(ID)
report = {'id': ID, 'files': {}, 'roundtrip': {}}
sd = f'{TREE}/sourcedata/{SLUG}'
sums = []
for n, exp in sorted(man.items()):
    dg = copy_verified(f'{SRC}/{n}', f'{sd}/{n}', exp)
    report['files'][f'sourcedata/{SLUG}/{n}'] = dg
    sums.append((n, dg))
write_sha256sums(sd, sums)
log('sourcedata done')

m = sio.loadmat(f'{SRC}/sEEG_data.mat')
freqs = m['frequencies']
fl = [float(x) for x in freqs.ravel()]
po = m['parcel_order']
parcel = {}
for row in po:
    code = int(np.asarray(row[0]).ravel()[0]); name = str(np.asarray(row[1]).ravel()[0])
    parcel[code] = name
os.makedirs(f'{TREE}/group/seeg', exist_ok=True)
with open(f'{TREE}/group/seeg/parcel_order.tsv', 'w') as f:
    f.write('parcel_code\tparcel_name\n')
    for row in po:
        f.write(f'{int(np.asarray(row[0]).ravel()[0])}\t{str(np.asarray(row[1]).ravel()[0])}\n')
with open(f'{TREE}/group/seeg/frequencies.tsv', 'w') as f:
    f.write('frequency_index\tfrequency_hz\n')
    for i, x in enumerate(fl):
        f.write(f'{i + 1}\t{x:g}\n')

FMT = lambda x: repr(float(x))
pt = []
for i in range(68):
    s = f'{i + 1:02d}'
    d = f'{TREE}/sub-{s}/ieeg'
    os.makedirs(d, exist_ok=True)
    cells = {
        'PLV': m['PLV'][i, 0], 'DFA': m['DFA'][i, 0], 'elec_dist': m['elec_dist'][i, 0], 'ref_mask': m['ref_masks'][i, 0],
        'acceptcontHealthy': m['acceptcontHealthy'][i, 0], 'acceptcontEZ': m['acceptcontEZ'][i, 0],
        'parcel_assign': m['parcel_assign'][0, i], 'HealthyEZcriteria': m['HealthyEZcriteria'][i:i + 1, :], 'frequencies': freqs,
    }
    n = cells['DFA'].shape[1]
    for k in ('elec_dist', 'ref_mask'):
        assert cells[k].shape == (n, n), (s, k, cells[k].shape)
    assert cells['PLV'].shape == (24, n, n)
    for k in ('acceptcontHealthy', 'acceptcontEZ'):
        assert cells[k].shape == (n, 1), (s, k)
    pa_ok = cells['parcel_assign'].shape == (n, 1)
    if not pa_ok:  # released as is (sub-18: 152 codes for 148 contacts); not aligned in the contacts TSV
        report.setdefault('parcel_assign_mismatch', {})[s] = {'n_contacts': n, 'parcel_assign_len': int(cells['parcel_assign'].shape[0])}
    out = f'{d}/sub-{s}_task-rest_desc-sync_connectivity.mat'
    sio.savemat(out, cells, format='5', do_compression=True, oned_as='column')
    back = sio.loadmat(out)
    ok = all(cells[k].dtype == back[k].dtype and cells[k].shape == back[k].shape and np.array_equal(cells[k], back[k], equal_nan=True) for k in cells)
    assert ok, out
    report['roundtrip'][os.path.relpath(out, TREE)] = {'exact': True, 'n_contacts': n}
    # DFA TSV
    tsv = f'{d}/sub-{s}_task-rest_desc-dfa_contacts.tsv'
    with open(tsv, 'w') as f:
        f.write('contact_index\t' + '\t'.join(f'dfa_{x:g}Hz' for x in fl) + '\n')
        for c in range(n):
            f.write(f'{c + 1}\t' + '\t'.join('n/a' if np.isnan(v) else FMT(v) for v in cells['DFA'][:, c]) + '\n')
    rows = [l.rstrip('\n').split('\t') for l in open(tsv)][1:]
    arr = np.array([[np.nan if v == 'n/a' else float(v) for v in r[1:]] for r in rows]).T
    assert np.array_equal(arr, cells['DFA'], equal_nan=True), tsv
    # contacts TSV
    ctsv = f'{d}/sub-{s}_task-rest_contacts.tsv'
    with open(ctsv, 'w') as f:
        f.write('contact_index\tparcel_code\tparcel_name\tnon_ez_contact\tez_contact\n')
        for c in range(n):
            pc = int(cells['parcel_assign'][c, 0]) if pa_ok else 'n/a'
            f.write(f"{c + 1}\t{pc}\t{parcel.get(pc, 'n/a')}\t{int(cells['acceptcontHealthy'][c, 0])}\t{int(cells['acceptcontEZ'][c, 0])}\n")
    rows = [l.rstrip('\n').split('\t') for l in open(ctsv)][1:]
    if pa_ok:
        assert [int(r[1]) for r in rows] == cells['parcel_assign'].ravel().tolist()
    assert [int(r[3]) for r in rows] == cells['acceptcontHealthy'].ravel().tolist()
    assert [int(r[4]) for r in rows] == cells['acceptcontEZ'].ravel().tolist()
    report['roundtrip'][os.path.relpath(tsv, TREE)] = {'exact': True}
    report['roundtrip'][os.path.relpath(ctsv, TREE)] = {'exact': True, 'unknown_parcel_codes': sorted({int(x) for x in cells['parcel_assign'].ravel()} - set(parcel)), 'parcel_assign_aligned': pa_ok}
    for p in (out, tsv, ctsv):
        report['files'][os.path.relpath(p, TREE)] = sha256(p)
    h = cells['HealthyEZcriteria'].ravel()
    pt.append((f'sub-{s}', n, int(cells['acceptcontHealthy'].sum()), int(cells['acceptcontEZ'].sum()), int(h[0]), int(h[1])))
del m
with open(f'{TREE}/participants.tsv', 'w') as f:
    f.write('participant_id\tsource_index\tn_contacts\tn_non_ez_contacts\tn_ez_contacts\tnon_ez_analysis\tez_analysis\n')
    for r in pt:
        yn = lambda b: 'yes' if b else 'no'
        f.write(f'{r[0]}\t{int(r[0][4:])}\t{r[1]}\t{r[2]}\t{r[3]}\t{yn(r[4])}\t{yn(r[5])}\n')
log('per-subject done')

report['texts'] = install_texts(TEXT, TREE)
os.makedirs(f'{TREE}/code', exist_ok=True)
for f in ('f05x_build056.py', 'f05x_lib.py'):
    shutil.copyfile(f'{ROOT}/harness/{f}', f'{TREE}/code/{f}')
report['validator'] = validate(TREE, REP)
log('validator', report['validator'])
pr = privacy_scan(TREE, REP)
report['privacy'] = {'n_findings': pr['n_findings'], 'stats': pr['stats'], 'findings_head': pr['findings'][:30]}
log('privacy', pr['n_findings'])
report['listing_n'] = len(tree_listing(TREE))
json.dump(report, open(REP + '/build.json', 'w'), indent=1)
print('RESULT', json.dumps({'validator': report['validator'], 'privacy_n': pr['n_findings'], 'n_files': report['listing_n'],
                            'roundtrip_all_exact': all(v['exact'] for v in report['roundtrip'].values())}))
