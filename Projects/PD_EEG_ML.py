"""Two-feature logistic regression with participant-level cross-validation.

Permutation testing uses scikit-learn's documented permutation_test_score:
https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.permutation_test_score.html
No feature selection, model search or threshold tuning is performed.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score, confusion_matrix, ConfusionMatrixDisplay,
    roc_auc_score, roc_curve,
)
from sklearn.model_selection import StratifiedKFold, permutation_test_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_DIR / "machine_learning"
FEATURES = ["Relative theta power", "Dominant frequency 4-13 Hz"]
SEED = 42
N_PERMUTATIONS = 1000


def main():
    data = pd.read_csv(PROJECT_DIR / "EEG_features.csv")
    X = data[FEATURES].to_numpy()
    y = data["Group"].map({"Control": 0, "PD": 1}).to_numpy()
    OUTPUT_DIR.mkdir(exist_ok=True)

    # The scaler learns only from the training participants in each fold.
    # Balanced weights give the two groups equal total weight during fitting.
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=1.0, solver="lbfgs", class_weight="balanced",
                           max_iter=1000, random_state=SEED),
    )
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    probabilities = np.zeros(len(y))
    predictions = np.zeros(len(y), dtype=int)
    fold_numbers = np.zeros(len(y), dtype=int)
    fold_results = []

    for fold, (train, test) in enumerate(cv.split(X, y), start=1):
        fitted = clone(model).fit(X[train], y[train])
        probabilities[test] = fitted.predict_proba(X[test])[:, 1]
        predictions[test] = fitted.predict(X[test])  # Fixed default threshold.
        fold_numbers[test] = fold
        fold_results.append({
            "Fold": fold, "Training N": len(train), "Test N": len(test),
            "Balanced accuracy": balanced_accuracy_score(y[test], predictions[test]),
            "ROC AUC": roc_auc_score(y[test], probabilities[test]),
        })
        print(f"Completed fold {fold}/5", flush=True)

    print(f"Running {N_PERMUTATIONS} label permutations...", flush=True)
    # Pass the splitter, not fixed label-dependent splits: stratification is
    # recomputed for each shuffled label set. The full pipeline is refitted.
    score, null_scores, pvalue = permutation_test_score(
        model, X, y, scoring="balanced_accuracy", cv=cv,
        n_permutations=N_PERMUTATIONS, random_state=SEED, n_jobs=1,
    )
    folds = pd.DataFrame(fold_results)
    matrix = confusion_matrix(y, predictions, labels=[0, 1])
    tn, fp, fn, tp = matrix.ravel()
    metrics = {
        "Participants": len(y), "PD N": int(y.sum()), "Control N": int((y == 0).sum()),
        "Mean fold balanced accuracy (permutation statistic)": score,
        "Pooled held-out balanced accuracy": balanced_accuracy_score(y, predictions),
        "Mean fold ROC AUC": folds["ROC AUC"].mean(),
        "Pooled held-out ROC AUC": roc_auc_score(y, probabilities),
        "PD sensitivity": tp / (tp + fn), "Control specificity": tn / (tn + fp),
        "Permutation p-value": pvalue, "Permutations": N_PERMUTATIONS,
        "Null mean balanced accuracy": null_scores.mean(), "Seed": SEED,
    }
    pd.DataFrame([metrics]).to_csv(OUTPUT_DIR / "model_metrics.csv", index=False)
    folds.to_csv(OUTPUT_DIR / "fold_scores.csv", index=False)
    pd.DataFrame({
        "Subject": data["Subject"], "Actual group": data["Group"],
        "Fold": fold_numbers, "Held-out PD score": probabilities,
        "Predicted group": np.where(predictions == 1, "PD", "Control"),
    }).to_csv(OUTPUT_DIR / "held_out_predictions.csv", index=False)
    pd.DataFrame({"Mean fold balanced accuracy": null_scores}).to_csv(
        OUTPUT_DIR / "permutation_scores.csv", index=False,
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    ConfusionMatrixDisplay(matrix, display_labels=["Control", "PD"]).plot(
        ax=axes[0], cmap="Blues", colorbar=False,
    )
    axes[0].set_title("Held-out predictions")
    fpr, tpr, _ = roc_curve(y, probabilities)
    axes[1].plot(fpr, tpr, color="#24748C",
                 label=f"Pooled AUC = {roc_auc_score(y, probabilities):.3f}")
    axes[1].plot([0, 1], [0, 1], "--", color="gray")
    axes[1].set(xlabel="False-positive rate", ylabel="True-positive rate",
                title="Held-out ROC curve")
    axes[1].legend(loc="lower right")
    axes[2].hist(null_scores, bins=25, color="#8AB7C3", edgecolor="white")
    axes[2].axvline(score, color="#C5653C", linewidth=2,
                    label=f"Observed = {score:.3f}")
    axes[2].set(xlabel="Mean fold balanced accuracy", ylabel="Label shuffles",
                title=f"1,000 permutations · p = {pvalue:.4f}")
    axes[2].legend()
    fig.suptitle("Logistic regression: theta power + dominant frequency", fontsize=15)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "classification_results.png", dpi=200)
    fig.savefig(OUTPUT_DIR / "classification_results.svg")
    plt.close(fig)
    for name, value in metrics.items():
        print(f"{name}: {value}")
    print("Confusion matrix (rows actual, columns predicted; Control, PD):\n", matrix)
    print(f"Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
