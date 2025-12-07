import numpy as np


def metrics(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float]:
    yt = y_true
    yp = y_pred
    yt_center = yt - yt.mean()
    yp_center = yp - yp.mean()
    denom = np.sqrt((yt_center**2).sum() * (yp_center**2).sum()) + 1e-12
    pearson = float((yt_center @ yp_center) / denom)
    rmse = np.sqrt(np.mean((yt - yp) ** 2))
    relrmse = float(rmse / (np.mean(np.abs(yt)) + 1e-12))
    return pearson, relrmse
