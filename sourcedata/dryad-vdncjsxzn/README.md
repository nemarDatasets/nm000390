# Data from: Brain criticality predicts individual levels of inter-areal synchronization in human electrophysiological data

Minimal dataset comprising results of phase synchronization and detrended fluctuation analysis which were computed for MEG data recorded from 52 healthy participants and SEEG data recorded 68 drug-resistant focal epileptic patients, and supporting files, as well as simulated data.

## MEG_data contents

wPLI: wPLI connectivity matrices for [192 sessions * 24 frequencies * 400*400 parcels]
DFA: DFA Exponents for [192 sessions * 24 frequencies * 400 parcels] 
efidmask: edge masks based on edge fidelity & cpPLV for [192 sessions * 400*400 parcels] 
sfidmask: masks based on rejected edges per parcel for [192 sessions * 400 parcels]
parcels: 400 parcel names of Yeo/Schaefer 7 networks parcellation (original naming convention)
frequencies: frequencies in Hz

## sEEG_data contents

PLV: complex PLV connectivity matrices for [68 subjects * 24 frequencies * number of contacts * number of contacts]
DFA: DFA Exponents for [68 subjects * 24 frequencies * number of contacts] 
HealthyEZcriteria: 68 subjects * 2 logical selection masks (1st column non-EZ subjects, 2nd column EZ)]
acceptcontHealthy: masks of non-EZ contacts for 68 subjects
acceptcontEZ: masks of EZ contacts for 68 subjects
ref_masks: additional reference selection masks for 68 subjects
elec_dist: electrodes (contacts) distances for 68 subjects
parcel_assign: contact morphing operators (parcel assignment codes)
parcel_order: ordered Yeo/Schaefer 17 networks cortical parcellation and subcortical regions names (parcel assignment codes starting with 0-Unknown)
frequencies: frequencies in Hz

## Simulated data contents

realistic_individuals_dense_stats.pickle - Contains average amplitude, DFA and variance of amplitude for each channel and PLV for each pair of channels, obtained with simulations using structural connectome with log-transformation. 

realistic_individuals_dense_noisy_stats.pickle - Contains average amplitude, DFA and variance of amplitude for each channel and PLV for each pair of channels, obtained with simulations using structural connectome with log-transformation and additional noise. 

uniform_individuals_dense_stats.pickle - Contains average amplitude, DFA and variance of amplitude for each channel and PLV for each pair of channels, obtained with simulations using uniform connectome with log-transformation. 

uniform_individuals_dense_noisy_stats.pickle - Contains average amplitude, DFA and variance of amplitude for each channel and PLV for each pair of channels, obtained with simulations using uniform connectome with log-transformation and additional noise. 

reg_individual_stats.pickle - Contains average amplitude, DFA and variance of amplitude for each channel and PLV for each pair of channels used to compute DFA and PLV regressions. 

non_log_stats.pickle - Contains average amplitude, DFA and variance of amplitude for each channel and PLV for each pair of channels, obtained with simulations using structural connectome without log-transformation.  

## Methods:

### MEG data acquisition

We recorded MEG data from 52 healthy participants (age: 31 ± 9.2, 27 male) during a 10-minute eyes-open resting-state session with a Vectorview/Triux (Elekta-Neuromag/MEGIN, Helsinki, Finland) 306-channel system (204 planar gradiometers and 102 magnetometers) at the Bio-Mag Laboratory, HUS Medical Imaging Center, Helsinki. Overall, 192 sessions of MEG data were obtained, with participants contributing on average 3.7 ± 4 sessions each. Participants were instructed to focus on a cross on the center of the screen in front of them. Bipolar horizontal and vertical electrooculography (EOG) were recorded for the detection of ocular artifacts. MEG and EOG were recorded at 1 kHz sampling rate. For each participant, T1-weighted anatomical MRI scans (MP-RAGE) at a resolution of 1 x 1 x 1 mm were obtained at Helsinki University Central Hospital with a 1.5-Tesla MRI scanner (Siemens, Munich, Germany) for head models and cortical surface reconstruction. The study protocol for MEG and MRI data was approved by the Coordinating Ethical Committee of Helsinki University Central Hospital (ID 290/13/03/2013), written informed consent was obtained from each participant prior to the experiment, and all research was carried out according to the Declaration of Helsinki.

### Cortical parcellation and source model

Volumetric segmentation of MRI data, flattening, cortical parcellation, and neuroanatomical labeling with the 400-parcel Schaefer atlas (Schaefer et al., 2018) was carried out with the FreeSurfer software (Fischl, 2012). The MNE software (Gramfort et al., 2014) was then used to create cortically constrained source models, for MEG-MRI colocalization, and for the preparation of the forward and inverse operators. The source models had dipole orientations fixed to pial-surface normals and a 5-mm interdipole separation throughout the cortex, which yielded around 5000-8000 source vertices per hemisphere.

### MEG data preprocessing and filtering

Temporal signal space separation (tSSS)  (Taula & Simola J, 2006) in the Maxfilter software (Elekta-Neuromag) was used to suppress extracranial noise from MEG sensors and to interpolate bad channels. We used independent components analysis adapted from the Fieldtrip toolbox (Oostenveld et al., 2011) to extract and identify components that were correlated with ocular artifacts (identified using the EOG signal), heartbeat artifacts (identified using the magnetometer signal as a reference), or muscle artifacts.

### MEG source localization

MNE software was used to create cortically constrained source models with 5-mm inter-dipole separation, for MEG-MRI colocalization, and for the preparation of the forward and inverse operators. We computed noise covariance matrices (NCMs) using the preprocessed MEG data filtered with finite-impulse-response (FIR) filters at 151-249 Hz, averaged across 10 s time-windows. NCMs were used for creating one inverse operator per session and the dSPM method with regularization parameter Î = 0.11. We then estimated vertex fidelity to obtain fidelity-weighted inverse operators that reduce the effects of spurious connections resulting from source leakage, and collapsed the inverse-transformed source time series into parcel time series in a manner that maximizes the source-reconstruction accuracy as in (Korhonen et al., 2014; Rouhinen et al., 2020; Siebenhühner et al., 2020). For each parcel pair (edge) we also computed cross-parcel phase-locking of the reconstructed simulated time-series, reflecting cross-parcel signal mixing, and excluded parcels and edges with low fidelity and high cross-parcel phase-locking, using individual thresholds to retain for each subject the top 90% parcels by fidelity and the bottom 95% of edges by cross-parcel mixing (14.9 ± 0.2% of parcels and 14.1 ± 0.1% of edges rejected on average per set).

### Acquisition of SEEG Data

We recorded ten minutes of stereo-EEG (SEEG) neuronal signals from 68 drug-resistant focal epileptic patients (age: 30 ± 9.4, 38 male) during the clinical assessment of the epileptogenic zone (EZ) at the Claudio Munari Epilepsy Surgery Centre in the Niguarda Ca' Granda Hospital, Milan. These resting-state recordings were free of seizure activity, and there were no seizures within one hour prior to or after the recording. During the recording, the patient was asked to lay down in their bed in a quiet resting-state with their eyes closed. Visual inspection of delta amplitude profiles from at least one EEG pair (C3-P3) and video recordings of the patient was used by trained technicians to ascertain the absence of signs of sleep or drowsiness. Intracranial monopolar (with contacts sharing the reference to a single white-matter contact) local-field potentials were acquired from brain tissue with platinum-iridium multi-lead electrodes. Between 8 to 15 contacts, each 2 mm long, 0.8 mm thick and with an inter-contact border-to-border distance of 1.5 mm (DIXI medical, Besancon, France), were present in each penetrating shaft, with the amounts of electrodes and their anatomical positions varying according to surgical requirements (Cardinale et al., 2013). 

Each subject had 17 ± 3 (range 9-23) shafts with a total of 153 ± 20 electrode contacts on average. The electrode positions were localized after implantation using CT scans and the SEEGA automatic contact localization. Structural MRIs were recorded before implantation and colocalized with post-implant CT scans using rigid-body coregistration (Arnulfo, Narizzano, et al., 2015). Individual patients' contacts were assigned to parcels of the Schaefer atlas (Schaefer et al., 2018). We acquired an average of 10 min of uninterrupted spontaneous resting-state activity with eyes closed with a 192-channel SEEG amplifier system (Nihon-Kohden Neurofax EEG-1100) at a sampling rate of 1 kHz. All patients were taking antiseizure medications (antiepileptic drugs (AEDs) with a large variation in the dosage and compounds, and the time elapsed the last drug administration and the SEEG resting-state recording was not controlled. Patients gave written informed consent for participation in research studies and for publication of results pertaining to their data. The ethical committee of the Niguarda Hospital, Milan, approved this study (ID 939) which was performed according to the Declaration of Helsinki.

### Filtering and preprocessing of SEEG data

The EZ is clinically defined as the brain regions where ictal activity initiates and propagates (Fisher et al 2014). In this work, the EZ in individual patient brains were identified and stringently confirmed by clinicians using peri-ictal and ictal SEEG recordings (Cossu et al., 2015). Analyses in Figures 2-3 included only contacts from tentatively healthy regions (nEZ), i.e., those in which ictal activity was not observed during SEEG monitoring, while analysis in Figure 4 were performed for EZ contacts. Subjects who had more than half of their contacts in the EZ (8 subjects) were excluded from the analysis of non-EZ data and included to EZ data. For non-EZ data, we included data from 57 subjects, with a total of 4453 non-EZ contacts (average per subject 78 ± 19, range 41-124). For EZ, we used only the contacts defined as being epileptogenic, i.e., contacts that were located within the EZ or were part of the seizure propagation network and all contacts from 8 subjects in whom > 50 % contacts had previously been classified as EZ. Subjects who had < 11 EZ contacts were excluded. Thus, both analyses had 57 subjects, with 49 common to both non-EZ and EZ analyses, 8 included only in non-EZ, and 8 only included in EZ analysis. The total number of EZ contacts was 1725 (average per subject 30 ± 17.5, range 11-79). The average distance between EZ and non-EZ contacts was similar. 

Defective contacts that demonstrated non-physiological activity (1.3 ± 1.2, range 0-10) were excluded from analyses, and 3 subjects in which more than 50% of contacts were defective were discarded from further analyses. As occasional inter-ictal events characterized as large amplitude spikes or sharp waves with wide spectral and spatial spread may bias the LRTC and phase synchronization estimates (Monto et al., 2007), we rejected temporal segments with such activity. For a full description see (Arnulfo et al., 2002). Briefly, we partitioned the signals in multiple non-overlapping time windows of 500ms and decomposed the full stack of time signals for each channel in 24 frequency bins. With such time-frequency decomposition, we first computed the mean and the standard deviation of the signal amplitude for each individual cortical channel and for each frequency. A single 500ms time window was considered containing putatively epileptic (i.e. high amplitude artifacts) if 10% of cortical contacts showed amplitude bursts stronger than their mean amplitude plus 5 times the standard deviation, and if the burst frequency spectrum involved at least 25% of the analyzed frequency bins.

We then referenced SEEG electrodes in grey matter to the closest contacts in white matter which yields signals with more accurate phase estimates (Arnulfo et al., 2015), FIR-filtered broad-band contact time series with a cutoff at 440 Hz, and removed 50 Hz line noise and its harmonics with notch FIR filters (Arnulfo et al., 2020).

### Analysis of inter-areal synchronization

Both cleaned MEG and SEEG data were filtered into complex-valued narrowband time series using Morlet wavelets (m = 5 width) with 24 logarithmically increasing center frequencies ranging from 2 to 225 Hz. We computed pairwise phase synchronization for all narrow-band frequencies between all 400 cortical parcels in source-reconstructed MEG data and between all non-epileptogenic contact-pairs of SEEG data for all narrow-band frequencies, using both the Phase-Locking Value (PLV) (Lachaux et al., 1999) and the weighted Phase Lag Index (wPLI) (Vinck et al., 2011) which, unlike PLV, is insensitive to zero-phase lagged synchronization caused by source-mixing in MEG.

### Detrended Fluctuation Analysis  

We used Detrended Fluctuation Analysis (DFA) (Hardstone et al., 2012; Linkenkaer-Hansen et al., 2001) to estimate monofractal scaling exponents of neuronal LRTCs that typically vary between 0.5 and 1. DFA was carried out in the Fourier domain (Nolte et al., 2019) with a Gaussian weight function used for detrending and using 25 log-linear windows from 5 s - 56 s. The fluctuations were fitted with a robust linear regression with a bisquare weight function to obtain the DFA exponents, all with negligible fit error. DFA exponents were computed for all contacts of all SEEG subjects and all parcels of all MEG sets, and for each narrow-band frequency.

### Data simulation with a hierarchical Kuramoto model:

We adapted a two-layer nested Kuramoto model of local and large-scale neuronal oscillatory dynamics with a hierarchical variant (Siebenhühner et al., 2020; Simola et al., 2022) consisting of 400 nodes (corresponding to 400 brain areas), each containing a conventional Kuramoto population of oscillators, N = 500. Oscillators were modeled with a gaussian frequency distribution with mean value of 10 and standard deviation of 3 to obtain properties of an oscillating network with coupled interactions, where the oscillation frequency per se does not play a role. The phases of all oscillators in a node were averaged to derive a single time series per region, whose absolute values were then taken to obtain the Kuramoto order parameter, corresponding to the amplitude envelope of local oscillations. In this model, a local control parameter K determines the homogeneous coupling strength within the nodes while a global control parameter L determines the inter-regional coupling strength based on the pairwise-connectivity white-matter axonal fiber counts estimated with structural diffusion tensor imaging (DTI) and obtained from the Human Connectome Project (HCP). White noise was added to the final timeseries to simulate device and environmental noise. This model yields dynamics both at local and large-scale network levels that are comparable to the empirically observable dynamics within and between cortical regions. Thus, we were able to model both local interactions within smaller regions (corresponding to cortical MEG parcels or SEEG contacts), and synchronization (corresponding to that between parcels or contacts across the whole brain). We assessed LRTCs for each region using DFA, and synchronization between all pairs of regions using the PLV (see below). 
