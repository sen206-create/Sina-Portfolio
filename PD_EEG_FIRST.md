# Parkinson’s Disease EEG: My First EEG Preprocessing Pipeline

> WORK IN PROGRESS — last updated 2 October 2026**
>
> I’m building this project as part of my computational neuroscience learning journey. This page documents the current version of my EEG preprocessing and quality-control code, including the decisions I’m learning to make along the way. The pipeline is still being tested and improved; it is not a finished or validated machine-learning system.
>
> **Come back to follow my progress.**
> I’ll update this write-up as I improve the code, process more participants, and move towards feature extraction and machine learning.

[View the Python script](PD_EEG_FIRST.py)

## What I’m building

The aim is to turn resting-state EEG recordings into cleaned epochs and a clear quality-control summary for each participant. Before comparing groups or training a model, I need to understand what has been removed, what has been reconstructed, and which recordings deserve further review.

The current script, `PD_EEG_FIRST.py`, combines automated bad-channel detection, independent component analysis (ICA), automatic component labelling, epoch repair, and an Excel summary.

![](Figures/EEGfigure.png)


## The dataset

I’m working with Rest eyes open, available as [OpenNeuro ds004584](https://openneuro.org/datasets/ds004584) and its [NEMAR copy, on004584](https://nemar.org/dataset/on004584). The dataset description reports 149 participants: 100 with Parkinson’s disease and 49 controls, recorded using a 64-channel BrainVision system. The published protocol describes two minutes of eyes-open rest. My first recording has 63 channels and a sampling frequency of 500 Hz.

The recordings come from the dataset authors; my contribution here is developing and documenting the processing code. Group comparisons and disease classification have not yet been implemented.

## Tools used

| Tool | Role in this project |
| --- | --- |
| MNE-Python | Read EEG, reference and filter signals, fit ICA, interpolate channels, and create epochs |
| PyPREP / `NoisyChannels` | Detect potentially unreliable electrodes |
| MNE-ICALabel | Predict the type of activity represented by each ICA component |
| AutoReject | Detect problematic epochs and repair channel segments |
| NumPy | Calculate interpolation rates and other numerical summaries |
| pandas | Convert participant summaries into a table |
| openpyxl | Write and format the Excel workbook |
| pathlib | Handle output paths and create folders |


## How the pipeline works

```text
Load EEG
   ↓
Detect and mark bad channels
   ↓
Average reference
   ↓
Create 1–40 Hz analysis data and 1–100 Hz ICA data
   ↓
Select epochs for ICA using an initial AutoReject pass
   ↓
Fit extended Infomax ICA → classify components with ICLabel
   ↓
Remove selected components from the analysis data
   ↓
Interpolate initially bad channels
   ↓
Create analysis epochs → final AutoReject pass
   ↓
Return results and write the QC spreadsheet
```

### Current review rules

A recording is flagged when any of these conditions holds:

- More than 20% of channels are initially marked bad.
- More than 30% of fitted ICA components are removed.
- Any channel is interpolated in more than 50% of epochs in the final pass.
- Less than 60 seconds of epochs remain.
- ICA reaches its iteration limit.

These are provisional screening thresholds, not validated subject-exclusion criteria.
A flag prompts review; a lack of flags does not certify clean data.

## What works now, and what I’m building next

### Implemented in the current code

- Bad-channel detection and average referencing.
- Separate analysis and ICA filters.
- Epoch selection before extended Infomax ICA.
- ICLabel-based component selection using a confidence threshold.
- Channel interpolation and final AutoReject processing.
- Per-participant QC dictionaries and a formatted Excel export.
- A subject loop

### Still in progress

- Save cleaned epochs, ICA solutions, and rejection logs per participant.
- Record failed participants and continue processing the remaining files.
- Validate the screening thresholds across a varied sample of recordings.
- Make the interpolation-rate denominator and reconstructed-channel burden explicit.

## Where the project is heading

Once the preprocessing and QC workflow is dependable, I plan to extract spectral features such as theta, alpha, and beta power and investigate participant-level comparisons and machine learning.

**Come back for updates as I move from preprocessing one recording to a reproducible workflow across the dataset, and eventually to feature extraction and machine learning.**

## References

- [Dataset: Rest eyes open — NEMAR](https://nemar.org/dataset/on004584)
- [Dataset: OpenNeuro ds004584](https://openneuro.org/datasets/ds004584)
- [MNE: ICA documentation](https://mne.tools/stable/generated/mne.preprocessing.ICA.html)
- [MNE-ICALabel: automated component classification](https://mne.tools/mne-icalabel/stable/generated/examples/00_iclabel.html)
- [AutoReject: preprocessing workflow with ICA](https://autoreject.github.io/stable/auto_examples/plot_autoreject_workflow.html)
