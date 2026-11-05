from pathlib import Path
import json

import numpy as np
import pandas as pd

from sklearn.base import clone
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier


PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_DIR / "machine_learning_model_comparison"

SEED = 42

FEATURE_SETS = {
    "theta_only": ["Relative theta power"],
    "both": [
        "Relative theta power",
        "Dominant frequency 4-13 Hz"
    ]
}

MODELS = {
    "logistic_regression": {
        "pipeline": Pipeline([
            ("scaler", StandardScaler()),
            ("classifier", LogisticRegression(
                solver="lbfgs",
                class_weight="balanced",
                max_iter=1000,
                random_state=SEED
            ))
        ]),
        "param_grid": {
            "classifier__C": [0.001, 0.01, 0.1, 1, 10, 100]
        },
        "score_type": "probability"
    },

    "rbf_svm": {
        "pipeline": Pipeline([
            ("scaler", StandardScaler()),
            ("classifier", SVC(
                kernel="rbf",
                class_weight="balanced"
            ))
        ]),
        "param_grid": {
            "classifier__C": [0.1, 1, 10],
            "classifier__gamma": ["scale", 0.1, 1]
        },
        "score_type": "decision_function"
    },

    "random_forest": {
        "pipeline": Pipeline([
            ("classifier", RandomForestClassifier(
                n_estimators=100,
                class_weight="balanced",
                random_state=SEED,
                n_jobs=1
            ))
        ]),
        "param_grid": {
            "classifier__max_depth": [2, 4, None],
            "classifier__min_samples_leaf": [2, 5]
        },
        "score_type": "probability"
    }
}


def main():
    data = pd.read_csv(PROJECT_DIR / "EEG_features.csv")

    labels = data["Group"].map({"Control": 0, "PD": 1})

    y = labels.to_numpy(dtype=int)

    OUTPUT_DIR.mkdir(exist_ok=True)

    cv_outer = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=SEED
    )

    outer_splits = list(cv_outer.split(data, y))
    comparison_results = []

    for feature_name, features in FEATURE_SETS.items():
        X = data[features].to_numpy(dtype=float)

        for model_name, config in MODELS.items():
            print(f"\n=== {feature_name} | {model_name} ===")

            combination_dir = OUTPUT_DIR / feature_name / model_name
            combination_dir.mkdir(parents=True, exist_ok=True)

            fold_results = []
            prediction_results = []

            for fold, (train_ix, test_ix) in enumerate(
                outer_splits, start=1
            ):
                X_train, X_test = X[train_ix], X[test_ix]
                y_train, y_test = y[train_ix], y[test_ix]

                cv_inner = StratifiedKFold(
                    n_splits=5,
                    shuffle=True,
                    random_state=SEED
                )

                grid_search = GridSearchCV(
                    estimator=clone(config["pipeline"]),
                    param_grid=config["param_grid"],
                    cv=cv_inner,
                    scoring="balanced_accuracy",
                    refit=True,
                    n_jobs=-1,
                    error_score="raise"
                )

                grid_search.fit(X_train, y_train)

                y_pred = grid_search.predict(X_test)

                if config["score_type"] == "decision_function":
                    y_score = grid_search.decision_function(X_test)
                else:
                    y_score = grid_search.predict_proba(X_test)[:, 1]

                outer_score = balanced_accuracy_score(
                    y_test, y_pred
                )
                outer_auc = roc_auc_score(y_test, y_score)

                fold_results.append({
                    "Feature set": feature_name,
                    "Model": model_name,
                    "Fold": fold,
                    "Best parameters": json.dumps(
                        grid_search.best_params_
                    ),
                    "Best inner balanced accuracy":
                        grid_search.best_score_,
                    "Outer balanced accuracy": outer_score,
                    "Outer ROC AUC": outer_auc
                })

                prediction_results.append(pd.DataFrame({
                    "Feature set": feature_name,
                    "Model": model_name,
                    "Subject": data.iloc[test_ix]["Subject"].to_numpy(),
                    "Fold": fold,
                    "Actual label": y_test,
                    "Predicted label": y_pred,
                    "PD score": y_score,
                    "Score type": config["score_type"]
                }))

                print(
                    f"Fold {fold}: "
                    f"balanced accuracy={outer_score:.3f}, "
                    f"AUC={outer_auc:.3f}, "
                    f"parameters={grid_search.best_params_}"
                )

            fold_table = pd.DataFrame(fold_results)
            prediction_table = pd.concat(
                prediction_results,
                ignore_index=True
            )

            fold_table.to_csv(
                combination_dir / "nested_fold_results.csv",
                index=False
            )

            prediction_table.to_csv(
                combination_dir / "nested_held_out_predictions.csv",
                index=False
            )

            comparison_results.append({
                "Feature set": feature_name,
                "Model": model_name,
                "Mean outer balanced accuracy":
                    fold_table["Outer balanced accuracy"].mean(),
                "SD outer balanced accuracy":
                    fold_table["Outer balanced accuracy"].std(ddof=1),
                "Mean outer ROC AUC":
                    fold_table["Outer ROC AUC"].mean()
            })

    comparison_table = pd.DataFrame(comparison_results)

    comparison_table.to_csv(
        OUTPUT_DIR / "model_comparison.csv",
        index=False
    )

    print("\n=== Model comparison ===")
    print(comparison_table.to_string(index=False))


if __name__ == "__main__":
    main()
