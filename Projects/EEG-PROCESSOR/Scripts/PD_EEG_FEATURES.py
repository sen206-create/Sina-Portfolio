import csv
from collections import Counter
from pathlib import Path

import mne
import numpy as np


PROJECT_DIR = Path(__file__).resolve().parent
PARTICIPANTS_PATH = PROJECT_DIR / "participants.tsv"
EPOCHS_DIR = PROJECT_DIR / "preprocessed_eeg"
OUTPUT_PATH = PROJECT_DIR / "EEG_features.csv"


def extract_features(epochs):
    picks = mne.pick_types(epochs.info, eeg=True, exclude="bads")
    n_samples = len(epochs.times)
    psd, freqs = epochs.compute_psd(
        method="welch", fmin=1.0, fmax=30.0, picks=picks,
        n_fft=n_samples, n_per_seg=n_samples, n_overlap=0,
        window="hamming", verbose=False,
    ).get_data(return_freqs=True)

    # Average epochs first, then calculate relative power per channel.
    channel_psd = psd.mean(axis=0)
    total_mask = (freqs >= 1) & (freqs < 30)
    theta_mask = (freqs >= 4) & (freqs < 8)
    peak_mask = (freqs >= 4) & (freqs <= 13)
    total_power = channel_psd[:, total_mask].sum(axis=1)
    if not np.isfinite(psd).all() or np.any(total_power <= 0):
        raise ValueError("Invalid spectral power.")
    # Equal-width frequency bins: the bin width cancels in the power ratio.
    theta = channel_psd[:, theta_mask].sum(axis=1) / total_power
    mean_psd = channel_psd.mean(axis=0)
    peak_index = int(np.argmax(mean_psd[peak_mask]))
    peak_freqs = freqs[peak_mask]
    return {
        "Relative theta power": float(theta.mean()),
        "Dominant frequency 4-13 Hz": float(peak_freqs[peak_index]),
        "Peak at search boundary": peak_index in (0, len(peak_freqs) - 1),
        "Clean epochs": len(epochs),
        "EEG channels used": len(picks),
        "Retained seconds": len(epochs) * n_samples / epochs.info["sfreq"],
    }


def main():
    with PARTICIPANTS_PATH.open(newline="") as stream:
        participants = list(csv.DictReader(stream, delimiter="\t"))

    rows = []
    for index, participant in enumerate(participants, start=1):
        subject = participant["participant_id"]
        print(f"[{index}/{len(participants)}] {subject}", flush=True)
        epochs = mne.read_epochs(
            EPOCHS_DIR / f"{subject}-epo.fif", preload=False, verbose=False,
        )
        features = extract_features(epochs)
        rows.append({"Subject": subject, "Group": participant["GROUP"], **features})

    with OUTPUT_PATH.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} participants: {dict(Counter(row['Group'] for row in rows))}")
    print(f"Output: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
