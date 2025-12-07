from pathlib import Path
import json
from typing import Dict, Iterable, List, Tuple

import numpy as np

from .fields import generate_1overf_field
from .features import build_feature_sets, build_feature_sets_with_topology
from .metrics import metrics
from .model import fit_ols, predict, train_test_split


def run_seed(
    seed: int,
    size: int,
    alpha: float,
    anisotropy: Tuple[float, float, float],
    patch_size: int,
    bins: int,
    models: Iterable[str],
    topo_thresholds: Iterable[float],
) -> Dict[str, Dict]:
    rng = np.random.default_rng(seed)
    rho = generate_1overf_field(size=size, alpha=alpha, anisotropy=anisotropy, seed=seed)
    if any(name.startswith("topo_") for name in models):
        feature_sets, y = build_feature_sets_with_topology(rho, patch_size=patch_size, bins=bins, topo_thresholds=topo_thresholds)
    else:
        feature_sets, y = build_feature_sets(rho, patch_size=patch_size, bins=bins)
    out: Dict[str, Dict] = {}
    for name in models:
        X = feature_sets.get(name)
        if X is None:
            continue
        Xtr, Xte, ytr, yte = train_test_split(X, y, rng)
        coef = fit_ols(Xtr, ytr)
        pred = predict(Xte, coef)
        p, r = metrics(yte, pred)
        out[name] = {"pearson": p, "relrmse": r, "coef": coef.tolist()}
    return out


def aggregate(
    seeds: Iterable[int],
    size: int,
    alpha: float,
    anisotropy: Tuple[float, float, float],
    patch_size: int,
    bins: int,
    models: Iterable[str],
    topo_thresholds: Iterable[float],
):
    rows: List[Tuple[str, int, float, float, List[float]]] = []
    for s in seeds:
        res = run_seed(
            seed=s,
            size=size,
            alpha=alpha,
            anisotropy=anisotropy,
            patch_size=patch_size,
            bins=bins,
            models=models,
            topo_thresholds=topo_thresholds,
        )
        for k, v in res.items():
            rows.append((k, s, v["pearson"], v["relrmse"], v["coef"]))
    return rows


def summarize(rows):
    summary = {}
    for name in {r[0] for r in rows}:
        vals_p = [r[2] for r in rows if r[0] == name]
        vals_r = [r[3] for r in rows if r[0] == name]
        summary[name] = {
            "pearson_mean": float(np.mean(vals_p)),
            "pearson_std": float(np.std(vals_p)),
            "relrmse_mean": float(np.mean(vals_r)),
            "relrmse_std": float(np.std(vals_r)),
            "n": len(vals_p),
        }
    return summary


def run_experiment(
    seeds: Iterable[int] = (0, 1, 2),
    size: int = 32,
    alpha: float = 1.8,
    anisotropy: Tuple[float, float, float] = (1.0, 0.8, 0.6),
    patch_size: int = 3,
    bins: int = 16,
    models: Iterable[str] = ("baseline_m", "geom_m_grad2_lap", "ent_m_entropy_lap"),
    topo_thresholds: Iterable[float] = (0.35, 0.5, 0.65),
    output_path: Path | None = None,
):
    seeds = list(seeds)
    models = list(models)
    rows = aggregate(
        seeds=seeds,
        size=size,
        alpha=alpha,
        anisotropy=anisotropy,
        patch_size=patch_size,
        bins=bins,
        models=models,
        topo_thresholds=topo_thresholds,
    )
    summary = summarize(rows)
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "setup": {
                "size": size,
                "alpha": alpha,
                "anisotropy": anisotropy,
                "seeds": seeds,
                "patch_size": patch_size,
                "bins": bins,
                "models": models,
            },
            "runs": rows,
            "summary": summary,
        }
        output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return rows, summary


def format_summary(summary: Dict[str, Dict]) -> str:
    lines = []
    for name, stats in summary.items():
        lines.append(
            f"{name}: pearson={stats['pearson_mean']:.6f} ±{stats['pearson_std']:.6f} ; "
            f"relRMSE={stats['relrmse_mean']:.6f} ±{stats['relrmse_std']:.6f} (n={stats['n']})"
        )
    return "\n".join(lines)


def main() -> None:
    _, summary = run_experiment()
    print("# Entropie vs géométrie — 1/f anisotrope")
    print(format_summary(summary))
    out_path = Path(__file__).resolve().parents[2] / "outputs" / "entropy_results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    run_experiment(output_path=out_path)


if __name__ == "__main__":
    main()
