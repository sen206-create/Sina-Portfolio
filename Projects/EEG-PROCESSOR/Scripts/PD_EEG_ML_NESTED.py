from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, permutation_test_score, GridSearchCV
from sklearn.metrics import (
    balanced_accuracy_score, confusion_matrix, ConfusionMatrixDisplay,
    roc_auc_score, roc_curve,
)
from sklearn.linear_model import LogisticRegression
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.use("Agg")


PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_DIR / "machine_learning_nested"
FEATURES = ["Relative theta power", "Dominant frequency 4-13 Hz"]
SEED = 42


def main():
    data = pd.read_csv(PROJECT_DIR / "EEG_features.csv")
    X = data[FEATURES].to_numpy()
    y = data["Group"].map({"Control": 0, "PD": 1}).to_numpy()
    OUTPUT_DIR.mkdir(exist_ok=True)

    c_values = [0.001, 0.01, 0.1, 1, 10, 100]
    param_grid = {"logisticregression__C": c_values}

    cv_outer = StratifiedKFold(
        n_splits=5, shuffle=True, random_state=SEED
    )
    outer_scores = []
    fold_results = []
    prediction_results = []
    for fold, (train_ix, test_ix) in enumerate(cv_outer.split(X, y), start=1):
        X_train, X_test = X[train_ix, :], X[test_ix, :]
        y_train, y_test = y[train_ix], y[test_ix]
        cv_inner = StratifiedKFold(
            n_splits=5, shuffle=True, random_state=SEED
        )
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(
                solver="lbfgs", class_weight="balanced", max_iter=1000, random_state=SEED
            )
        )

        grid_search = GridSearchCV(
            model, param_grid, cv=cv_inner, scoring='balanced_accuracy', refit=True
        )
        grid_search.fit(X_train, y_train)

        y_pred = grid_search.predict(X_test)
        y_prob = grid_search.predict_proba(X_test)[:, 1]
        outer_score = balanced_accuracy_score(y_test, y_pred)
        outer_scores.append(outer_score)

        print(f"Fold {fold}")
        print("Best C:", grid_search.best_params_["logisticregression__C"])
        print("Best inner balanced accuracy:", grid_search.best_score_)

        print(f"Fold {fold}: outer balanced accuracy = {outer_score:.3f}")
        fold_results.append({
            "Fold": fold,
            "Best C": grid_search.best_params_["logisticregression__C"],
            "Best inner balanced accuracy": grid_search.best_score_,
            "Outer balanced accuracy": outer_score,
            "Outer ROC AUC": roc_auc_score(y_test, y_prob)
        })
        prediction_results.append(pd.DataFrame({
            "Subject": data.iloc[test_ix]["Subject"].to_numpy(),
            "Fold": fold,
            "Actual label": y_test,
            "Predicted label": y_pred,
            "PD score": y_prob
        }))
    print("Outer scores:", outer_scores)
    print(f"Mean outer balanced accuracy: {np.mean(outer_scores):.3f}")
    fold_table = pd.DataFrame(fold_results)

    prediction_table = pd.concat(
        prediction_results,
        ignore_index=True
    )

    fold_table.to_csv(
        OUTPUT_DIR / "nested_fold_results.csv",
        index=False
    )

    prediction_table.to_csv(
        OUTPUT_DIR / "nested_held_out_predictions.csv",
        index=False
    )


if __name__ == "__main__":
    main()
