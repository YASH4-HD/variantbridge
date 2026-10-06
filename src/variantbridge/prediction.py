"""EvoResist-AI exploratory prediction baseline (moved unchanged from v0.1 app.py).

Held-out Random Forest + permutation importance. Methodology intentionally NOT expanded in
v0.2. Predictive importance is not mechanism or causation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, classification_report, roc_auc_score)
from sklearn.model_selection import train_test_split


def holdout_random_forest(X: pd.DataFrame, y: pd.Series) -> dict:
    X = X.fillna(X.mean()).fillna(0)  # v0.1 behaviour preserved for the AMR prediction path only
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=.25, stratify=y, random_state=42)
    model = RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=42, n_jobs=-1).fit(Xtr, ytr)
    model.set_params(n_jobs=1)  # prediction only: avoids joblib thread overhead in the 140x6 importance re-predictions; results unchanged
    pr = model.predict(Xte)
    pb = model.predict_proba(Xte)[:, 1]
    pi = permutation_importance(model, Xte, yte, n_repeats=6, random_state=42, scoring="roc_auc")
    fi = pd.DataFrame({"feature": X.columns, "importance": pi.importances_mean, "sd": pi.importances_std}
                      ).sort_values("importance", ascending=False)
    return {"roc_auc": float(roc_auc_score(yte, pb)), "accuracy": float(accuracy_score(yte, pr)),
            "balanced_accuracy": float(balanced_accuracy_score(yte, pr)),
            "importance": fi, "report": pd.DataFrame(classification_report(yte, pr, output_dict=True)).T,
            "n_test": int(len(yte))}
