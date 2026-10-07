# Brain criticality predicts individual levels of inter-areal synchronization in human electrophysiological data (Fuscà, Siebenhühner, Wang et al., 2023): SEEG/MEG synchronization and DFA (derivative)

**This is a processed-data (derivative) dataset.** It repackages the authors' public Dryad release
(doi:[10.5061/dryad.vdncjsxzn](https://doi.org/10.5061/dryad.vdncjsxzn), version 7, CC0-1.0), described by the
authors as a "Minimal dataset comprising results of phase synchronization and detrended fluctuation analysis which were
computed for MEG data recorded from 52 healthy participants and SEEG data recorded 68 drug-resistant focal epileptic
patients, and supporting files, as well as simulated data." **The raw recordings are not part of the release and are
not included here.** The article's Data availability statement: "Raw electrophysiological data cannot be shared
publicly due to regulations imposed by the Ethical Committees but can be shared for collaborative efforts upon
request." (We found no open release of the raw 1 kHz resting SEEG of this cohort on Zenodo, OpenNeuro, DANDI or NEMAR.)

Article: Fuscà M, Siebenhühner F, Wang SH, Myrov V, Arnulfo G, Nobili L, Palva JM, Palva S (2023). Brain criticality
predicts individual levels of inter-areal synchronization in human electrophysiological data. *Nat Commun* 14:4736.
doi:[10.1038/s41467-023-40056-9](https://doi.org/10.1038/s41467-023-40056-9). Code: <https://github.com/palvalab/DFA_Synch/>.

## Data (article, Methods)

- **SEEG, 68 patients** (age 30 ± 9.4 years, 38 male; per-patient values not released): 10 min of eyes-closed resting
  state at the "Claudio Munari" Epilepsy Surgery Centre, Niguarda Ca' Granda Hospital, Milan, without seizures within
  one hour before or after. Monopolar recordings (shared white-matter reference), 8-15 contacts per shaft, 17 ± 3
  shafts and 153 ± 20 contacts per patient on average; contacts localised from CT (SEEGA). 500-ms windows with
  epileptic artifacts rejected; grey-matter contacts referenced to the closest white-matter contact; FIR broadband
  filtering (cutoff 440 Hz) and 50 Hz notch filters; 24 Morlet wavelets (m = 5), 2-225 Hz; phase-locking value (PLV)
  between contact pairs and DFA exponents of the amplitude envelopes per contact. Epileptogenic-zone (EZ) contacts were
  identified by clinicians from peri-ictal and ictal SEEG.
- **MEG, 52 healthy subjects** (age 31 ± 9.2 years, 27 male), **192 eyes-open resting-state sessions** (3.7 ± 4 per
  subject), Vectorview/Triux 306-channel, BioMag Laboratory, Helsinki; tSSS, ICA; source-reconstructed to the
  400-parcel Schaefer atlas; wPLI and DFA. The release gives no session-to-subject mapping, so the MEG data stay
  group-level.
- **Simulations**: hierarchical Kuramoto models (uniform vs realistic structural connectomes), pickled statistics.

## Layout

```
sub-NN/ieeg/                                   NN = 01..68 = row of the patient in sEEG_data.mat
  sub-NN_task-rest_desc-sync_connectivity.mat   that patient's arrays (MATLAB v5, values/dtypes unchanged):
      PLV [24 freq x C x C] complex64   complex phase-locking value between contacts
      DFA [24 x C]                      DFA exponents per contact
      elec_dist [C x C]                 contact-to-contact distances
      ref_mask [C x C]                  reference selection mask (release variable ref_masks)
      acceptcontHealthy [C x 1]         non-EZ contact mask
      acceptcontEZ [C x 1]              EZ contact mask
      parcel_assign [C x 1]             parcel code per contact (see group/seeg/parcel_order.tsv)
      HealthyEZcriteria [1 x 2]         patient selected for non-EZ (col 1) / EZ (col 2) analyses
      frequencies [1 x 24]              Hz
  sub-NN_task-rest_desc-dfa_contacts.tsv        DFA exponents, one row per contact, one column per frequency
  sub-NN_task-rest_contacts.tsv                 contact index, parcel code and name, non-EZ / EZ mask
group/seeg/parcel_order.tsv                    parcel codes and names: Yeo/Schaefer 17-network cortical parcels and
                                               subcortical regions (code 0 = Unknown), from parcel_order
group/seeg/frequencies.tsv                     the 24 frequencies (Hz)
sourcedata/dryad-vdncjsxzn/                    the 4 original files byte-for-byte (SHA256SUMS), incl. the authors' README.md
    sEEG_data.mat        all 68 patients in one file (cell arrays)
    MEG_data.mat         MATLAB v7.3 (HDF5): wPLI [400 x 400 x 24 x 192 sessions], DFA [400 x 24 x 192],
                         efidmask [400 x 400 x 192] edge masks (edge fidelity and cpPLV), sfidmask [400 x 192] parcel
                         masks, parcels (400 Schaefer 7-network names), frequencies (Hz)
    Simulated_data.zip   six pickles: average amplitude, DFA and amplitude variance per channel and PLV per channel pair
                         for the model variants (see the authors' README.md)
code/                                          f05x_build056.py, f05x_lib.py
```

Every per-subject derivative file has a JSON sidecar of the same name (description, `Sources` pointing to the original release file, variables, how it was repackaged).

`C` is the number of contacts of the patient (70-148). Exception as released: for `sub-18`, `parcel_assign` holds 152
codes for 148 contacts, so it cannot be aligned to the contacts; it is kept unchanged in the MAT file and the
`parcel_code`/`parcel_name` columns of that patient's `_contacts.tsv` are `n/a`. Some contacts carry parcel code `-1`,
which is not listed in `parcel_order`; their `parcel_name` is `n/a`. Contact order is the same in every array and table of a patient;
contact names and coordinates are not part of the release. The `.mat` files and the `_contacts.tsv` tables are not BIDS
raw data types and are listed in `.bidsignore`. `participants.tsv` gives, per patient, the number of contacts, the
size of the non-EZ and EZ masks and the two `HealthyEZcriteria` selection flags (40 patients selected for non-EZ, 42
for EZ analyses).

MEG_data.mat and Simulated_data.zip are group-level files without a per-subject structure; use them directly from
`sourcedata/dryad-vdncjsxzn/` (`h5py` or MATLAB for the v7.3 MAT file; the pickles need Python and should only be
unpickled from this trusted source).

Round-trip: every per-patient array was re-read and compared exactly (dtype, shape, values) with the original cell;
the TSV tables re-parse to the original values (`code/f05x_build056.py`).

## Loading

```python
import scipy.io as sio
m = sio.loadmat('sub-01/ieeg/sub-01_task-rest_desc-sync_connectivity.mat')
plv = abs(m['PLV'])          # 24 x C x C
nonez = m['acceptcontHealthy'].ravel().astype(bool)
```

## Privacy

The release identifies patients only by row index. Text files, MAT strings (parcel names) and zip member names were
scanned for names, dates, record numbers and paths; nothing identifying was found.

## Ethics approval

SEEG: "Patients gave written informed consent for participation in research studies and for publication of results
pertaining to their data. The ethical committee of the Niguarda Hospital, Milan, approved this study (ID 939) which was
performed according to the Declaration of Helsinki." (Fuscà et al. 2023, *Nat Commun* 14:4736, Methods, Acquisition of
SEEG data.)

MEG: "The study protocol for MEG and MRI data was approved by the Coordinating Ethical Committee of Helsinki University
Central Hospital (ID 290/13/03/2013), written informed consent was obtained from each subject prior to the experiment,
and all research was carried out according to the Declaration of Helsinki." (same article, Methods, Acquisition of MEG
and MRI data.)

## Funding

"This work was supported by grants from the Academy of Finland (SA 1266745, 1296304 to J.M.P. and SA 325404 to S.P.),
from the Finnish Cultural Foundation to S.H.W. (postdoc fellowship 00220071), and from the Sigrid Jusélius Foundation to
S.P. and J.M.P." (article, Acknowledgements)

## Related

The SEEG cohort of Siebenhühner et al. (2020, *PLoS Biol*, doi:10.1371/journal.pbio.3000685; Dryad
doi:10.5061/dryad.0k86k80; 59 patients) comes from the same Niguarda programme; neither release maps subjects between
the two.

## License and citation

CC0-1.0 (Dryad record license). Please cite the article (doi:10.1038/s41467-023-40056-9) and the Dryad dataset
(doi:10.5061/dryad.vdncjsxzn).
