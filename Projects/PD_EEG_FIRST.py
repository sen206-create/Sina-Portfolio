from pathlib import Path

import mne
import numpy as np
import pandas as pd
from autoreject import AutoReject
from mne.preprocessing import ICA
from mne_icalabel import label_components
from pyprep import NoisyChannels


DATA_DIR = Path("/Users/sinaenayati/Desktop/v1")
OUTPUT_DIR = Path(__file__).resolve().parent / "preprocessed_eeg"


def missing_channel_positions(raw):
    eeg_picks = mne.pick_types(raw.info, eeg=True, exclude=[])
    if len(eeg_picks) == 0:
        raise ValueError("No EEG channels found.")
    missing = []
    for pick in eeg_picks:
        position = raw.info["chs"][pick]["loc"][:3]
        if not np.isfinite(position).all() or np.all(position == 0):
            missing.append(raw.ch_names[pick])
    return missing


def preprocess_eeg(raw):
    detector = NoisyChannels(
        raw.copy(),
        random_state=42
    )

    detector.find_all_bads()
    bads = detector.get_bads()

    raw_marked = raw.copy()
    raw_marked.info["bads"] = sorted(
        set(raw.info["bads"]) | set(bads)
    )

    raw_marked.set_eeg_reference("average", projection=False)

    raw_filtered = raw_marked.copy().filter(l_freq=1, h_freq=40.0)
    raw_for_ica = raw_marked.copy().filter(l_freq=1, h_freq=100.0)

    epochs_for_ica = mne.make_fixed_length_epochs(
        raw_for_ica,
        duration=2,
        overlap=0,
        preload=True,
    )
    ar_before_ica = AutoReject(random_state=42)
    ar_before_ica.fit(epochs_for_ica)
    ica_epoch_reject_log = ar_before_ica.get_reject_log(epochs_for_ica)

    epochs_for_ica_fit = epochs_for_ica[~ica_epoch_reject_log.bad_epochs]

    if len(epochs_for_ica_fit) == 0:
        raise ValueError("No epochs remaining for ICA.")
    ica = ICA(
        n_components=None,
        method='infomax',
        fit_params=dict(extended=True),
        random_state=42,
        max_iter='auto',
    )
    ica.fit(epochs_for_ica_fit, picks="eeg")
    ic_labels = label_components(epochs_for_ica_fit, ica, method="iclabel")
    labels = ic_labels["labels"]
    confidence = ic_labels["y_pred_proba"]
    exclude_idx = [
        idx
        for idx, (label, score) in enumerate(zip(labels, confidence))
        if label not in ["brain", "other"] and score >= 0.90
    ]
    ica.exclude = exclude_idx
    raw_clean = ica.apply(raw_filtered.copy())
    raw_clean.interpolate_bads(reset_bads=True)

    epochs = mne.make_fixed_length_epochs(
        raw_clean,
        duration=4,
        overlap=0,
        preload=True
    )
    ar = AutoReject(random_state=42)
    epochs_clean, reject_log = ar.fit_transform(epochs, return_log=True)
    epochs_clean.set_eeg_reference("average", projection=False)

    if len(epochs_clean) == 0:
        raise ValueError("No clean epochs remaining.")

    return epochs_clean


def main():
    participants = pd.read_csv(DATA_DIR / "participants.tsv", sep="\t")
    if participants["participant_id"].isna().any() or not participants["participant_id"].is_unique:
        raise ValueError("Participant IDs must be present and unique.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records = []
    retry_ids = {
        f"sub-{number:03d}"
        for number in range(1, 150)
    }
    for subject_id in participants["participant_id"]:
        if subject_id not in retry_ids:
            continue
        print(f"Processing {subject_id}...", flush=True)
        record = {"Subject": subject_id, "Status": "Failed", "Reason": "",
                  "Clean epochs": 0, "Retained seconds": 0, "File": ""}
        try:
            file_path = DATA_DIR / subject_id / "eeg" / \
                f"{subject_id}_task-Rest_eeg.set"
            raw = mne.io.read_raw_eeglab(file_path, preload=True)
            extra_channels = [
                name for name in raw.ch_names
                if name.strip().casefold() in {"resp", "x", "y", "z"}
            ]

            if extra_channels:
                print(f"{subject_id}: omitting channels {extra_channels}")
                raw.drop_channels(extra_channels)
            # Rebuild position metadata after removing the extra channels.
            # Preserve the remaining EEG coordinates and remove invalid leftover points.
            positions = raw.get_montage().get_positions()

            positions["ch_pos"] = {
                raw.ch_names[pick]: raw.info["chs"][pick]["loc"][:3].copy()
                for pick in mne.pick_types(raw.info, eeg=True, exclude=[])
            }

            for key in ("hsp", "hpi"):
                points = positions[key]
                if points is not None:
                    points = points[np.isfinite(points).all(axis=1)]
                    positions[key] = points if len(points) else None

            raw.set_montage(mne.channels.make_dig_montage(**positions))
            missing = missing_channel_positions(raw)
            if missing:
                record["Status"] = "Excluded"
                record["Reason"] = "Missing EEG channel positions: " + \
                    ", ".join(missing)
                print(f"Excluded {subject_id}: {record['Reason']}", flush=True)
            else:
                epochs_clean = preprocess_eeg(raw)
                output_path = OUTPUT_DIR / f"{subject_id}-epo.fif"
                epochs_clean.save(output_path, overwrite=True)
                record.update({
                    "Status": "Saved",
                    "Clean epochs": len(epochs_clean),
                    "Retained seconds": len(epochs_clean) * 4,
                    "File": str(output_path),
                })
                print(f"Saved {output_path}", flush=True)
        except Exception as error:
            record["Reason"] = f"{type(error).__name__}: {error}"
            print(f"Failed {subject_id}: {record['Reason']}", flush=True)

        records.append(record)
        pd.DataFrame(records).to_csv(
            OUTPUT_DIR / "processing_log.csv", index=False)

    counts = pd.DataFrame(records)["Status"].value_counts() if records else {}
    print(f"Finished: {counts.get('Saved', 0)} saved, "
          f"{counts.get('Excluded', 0)} excluded, {counts.get('Failed', 0)} failed.")
    print(f"Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
