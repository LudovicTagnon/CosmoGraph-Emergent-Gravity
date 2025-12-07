import numpy as np


def train_test_split(X: np.ndarray, y: np.ndarray, rng: np.random.Generator, test_ratio: float = 0.4):
    idx = np.arange(len(y))
    rng.shuffle(idx)
    t = int(len(y) * (1 - test_ratio))
    train, test = idx[:t], idx[t:]
    return X[train], X[test], y[train], y[test]


def fit_ols(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    Xb = np.c_[X, np.ones(len(X))]
    coef, *_ = np.linalg.lstsq(Xb, y, rcond=None)
    return coef


def predict(X: np.ndarray, coef: np.ndarray) -> np.ndarray:
    return np.c_[X, np.ones(len(X))] @ coef
