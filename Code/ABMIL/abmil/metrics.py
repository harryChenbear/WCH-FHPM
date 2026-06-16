import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, f1_score, confusion_matrix


def _safe_auc(y_true, prob):
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(roc_auc_score(y_true, prob))


def _safe_auprc(y_true, prob):
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(average_precision_score(y_true, prob))


def binary_metrics(y_true, prob, threshold=0.5):
    y_true = np.asarray(y_true).astype(int)
    prob = np.asarray(prob).astype(float)
    pred = (prob >= threshold).astype(int)

    out = {
        "threshold": threshold,
        "accuracy": float(accuracy_score(y_true, pred)),
        "f1": float(f1_score(y_true, pred, zero_division=0)),
        "auc": _safe_auc(y_true, prob),
        "auprc": _safe_auprc(y_true, prob),
    }

    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    out["sensitivity"] = float(tp / (tp + fn)) if (tp + fn) else np.nan
    out["specificity"] = float(tn / (tn + fp)) if (tn + fp) else np.nan
    out["ppv"] = float(tp / (tp + fp)) if (tp + fp) else np.nan
    out["npv"] = float(tn / (tn + fn)) if (tn + fn) else np.nan
    out["tn"] = int(tn)
    out["fp"] = int(fp)
    out["fn"] = int(fn)
    out["tp"] = int(tp)

    return out


def bootstrap_ci(y_true, prob, threshold=0.5, n_bootstrap=1000, seed=42):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true).astype(int)
    prob = np.asarray(prob).astype(float)
    n = len(y_true)

    records = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        yt = y_true[idx]
        pr = prob[idx]
        if len(np.unique(yt)) < 2:
            continue
        records.append(binary_metrics(yt, pr, threshold=threshold))

    keys = ["auc", "auprc", "accuracy", "f1", "sensitivity", "specificity", "ppv", "npv"]
    out = {}

    for k in keys:
        vals = np.array([r[k] for r in records], dtype=float)
        vals = vals[~np.isnan(vals)]
        out[f"{k}_mean"] = float(np.mean(vals)) if len(vals) else np.nan
        out[f"{k}_low"] = float(np.percentile(vals, 2.5)) if len(vals) else np.nan
        out[f"{k}_high"] = float(np.percentile(vals, 97.5)) if len(vals) else np.nan

    return out
