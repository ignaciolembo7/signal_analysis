from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

import numpy as np
from scipy.optimize import curve_fit, least_squares


ModelFn = Callable[..., np.ndarray]


@dataclass(frozen=True)
class ParametricFit:
    values: dict[str, float]
    errors: dict[str, float]
    yhat: np.ndarray
    rmse: float
    chi2: float
    r2: float
    method: str
    raw_result: Any | None = field(default=None, repr=False, compare=False)


@dataclass(frozen=True)
class CurveFitParameter:
    name: str
    value: float
    lower: float
    upper: float
    vary: bool = True


def _as_float_array(values: Sequence[float] | np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=float)


def rmse(y: np.ndarray, yhat: np.ndarray) -> float:
    y_arr = _as_float_array(y)
    yhat_arr = _as_float_array(yhat)
    return float(np.sqrt(np.mean((y_arr - yhat_arr) ** 2)))


def chi2(y: np.ndarray, yhat: np.ndarray) -> float:
    y_arr = _as_float_array(y)
    yhat_arr = _as_float_array(yhat)
    return float(np.sum((y_arr - yhat_arr) ** 2))


def r2_score(y: np.ndarray, yhat: np.ndarray) -> float:
    y_arr = _as_float_array(y)
    yhat_arr = _as_float_array(yhat)
    ss_res = float(np.sum((y_arr - yhat_arr) ** 2))
    ss_tot = float(np.sum((y_arr - np.mean(y_arr)) ** 2))
    return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan


def rmse_log(y: np.ndarray, yhat: np.ndarray) -> float:
    y_safe = np.clip(_as_float_array(y), 1e-12, None)
    yhat_safe = np.clip(_as_float_array(yhat), 1e-12, None)
    return float(np.sqrt(np.mean((np.log(y_safe) - np.log(yhat_safe)) ** 2)))


def parameter_error_column(param_name: str) -> str:
    if param_name == "D0_m2_ms":
        return "D0_err_m2_ms"
    if param_name == "D0_mm2_s":
        return "D0_err_mm2_s"
    if param_name == "tc_ms":
        return "tc_err_ms"
    if param_name == "g0_mTm":
        return "g0_err_mTm"
    return f"{param_name}_err"


def _covariance_errors(param_names: Sequence[str], pcov: np.ndarray | None) -> dict[str, float]:
    if pcov is None or not np.all(np.isfinite(pcov)):
        return {parameter_error_column(name): np.nan for name in param_names}
    perr = np.sqrt(np.clip(np.diag(pcov), 0.0, np.inf))
    return {parameter_error_column(name): float(err) for name, err in zip(param_names, perr)}


def _least_squares_errors(
    *,
    param_names: Sequence[str],
    popt: np.ndarray,
    jacobian: np.ndarray,
    residual: np.ndarray,
    log_params: set[str],
) -> dict[str, float]:
    if jacobian.size == 0 or residual.size <= len(param_names):
        return {parameter_error_column(name): np.nan for name in param_names}

    try:
        jtj_inv = np.linalg.pinv(jacobian.T @ jacobian)
    except np.linalg.LinAlgError:
        return {parameter_error_column(name): np.nan for name in param_names}

    dof = max(1, int(residual.size) - int(len(param_names)))
    residual_var = float(np.sum(residual ** 2) / dof)
    cov_opt = jtj_inv * residual_var
    opt_err = np.sqrt(np.clip(np.diag(cov_opt), 0.0, np.inf))

    errors: dict[str, float] = {}
    for i, name in enumerate(param_names):
        err = float(opt_err[i])
        if not np.isfinite(err):
            errors[parameter_error_column(name)] = np.nan
        elif name in log_params:
            errors[parameter_error_column(name)] = float(abs(popt[i]) * err)
        else:
            errors[parameter_error_column(name)] = err
    return errors


def fit_curve_fit(
    model: ModelFn,
    x: np.ndarray,
    y: np.ndarray,
    *,
    param_names: Sequence[str],
    p0: Sequence[float],
    bounds: tuple[Sequence[float], Sequence[float]],
    sigma: np.ndarray | None = None,
    absolute_sigma: bool = False,
    maxfev: int = 200000,
    method: str = "scipy_curve_fit",
) -> ParametricFit:
    x_arr = _as_float_array(x)
    y_arr = _as_float_array(y)
    popt, pcov = curve_fit(
        model,
        x_arr,
        y_arr,
        p0=list(p0),
        bounds=(list(bounds[0]), list(bounds[1])),
        sigma=sigma,
        absolute_sigma=absolute_sigma,
        maxfev=int(maxfev),
    )
    popt = _as_float_array(popt)
    yhat = _as_float_array(model(x_arr, *popt))
    return ParametricFit(
        values={name: float(value) for name, value in zip(param_names, popt)},
        errors=_covariance_errors(param_names, pcov),
        yhat=yhat,
        rmse=rmse(y_arr, yhat),
        chi2=chi2(y_arr, yhat),
        r2=r2_score(y_arr, yhat),
        method=str(method),
    )


def fit_curve_fit_parameters(
    model_from_params: Callable[[dict[str, float]], np.ndarray],
    y: np.ndarray,
    *,
    parameters: Sequence[CurveFitParameter],
    x: np.ndarray | None = None,
    maxfev: int = 200000,
    method: str = "scipy_curve_fit",
    fixed_method: str = "fixed",
) -> ParametricFit:
    """Fit a curve by varying only the parameters marked as variable."""
    y_arr = _as_float_array(y)
    x_arr = np.zeros_like(y_arr) if x is None else _as_float_array(x)
    params = list(parameters)
    values = {p.name: float(p.value) for p in params}
    variable = [p for p in params if p.vary]

    if not variable:
        yhat = _as_float_array(model_from_params(values))
        return ParametricFit(
            values=values,
            errors={parameter_error_column(p.name): np.nan for p in params},
            yhat=yhat,
            rmse=rmse(y_arr, yhat),
            chi2=chi2(y_arr, yhat),
            r2=r2_score(y_arr, yhat),
            method=str(fixed_method),
        )

    def wrapped_model(x_model: np.ndarray, *fit_values: float) -> np.ndarray:
        del x_model
        local_values = dict(values)
        local_values.update({p.name: float(v) for p, v in zip(variable, fit_values)})
        return _as_float_array(model_from_params(local_values))

    fit = fit_curve_fit(
        wrapped_model,
        x_arr,
        y_arr,
        param_names=[p.name for p in variable],
        p0=[p.value for p in variable],
        bounds=([p.lower for p in variable], [p.upper for p in variable]),
        maxfev=maxfev,
        method=method,
    )
    values.update(fit.values)
    errors = {parameter_error_column(p.name): np.nan for p in params if not p.vary}
    errors.update(fit.errors)
    return ParametricFit(
        values=values,
        errors=errors,
        yhat=fit.yhat,
        rmse=fit.rmse,
        chi2=fit.chi2,
        r2=fit.r2,
        method=fit.method,
        raw_result=fit.raw_result,
    )


def fit_least_squares(
    model: ModelFn,
    x: np.ndarray,
    y: np.ndarray,
    *,
    param_names: Sequence[str],
    p0: Sequence[float],
    bounds: tuple[Sequence[float], Sequence[float]],
    log_params: Sequence[str] = (),
    sigma: np.ndarray | None = None,
    max_nfev: int = 200000,
    method: str = "least_squares",
) -> ParametricFit:
    x_arr = _as_float_array(x)
    y_arr = _as_float_array(y)
    names = [str(name) for name in param_names]
    log_set = {str(name) for name in log_params}

    opt_p0: list[float] = []
    opt_lower: list[float] = []
    opt_upper: list[float] = []
    for name, value, lower, upper in zip(names, p0, bounds[0], bounds[1]):
        if name in log_set:
            if float(value) <= 0 or float(lower) <= 0:
                raise ValueError(f"{name} must be positive when fitted in log-space.")
            opt_p0.append(float(np.log(value)))
            opt_lower.append(float(np.log(lower)))
            opt_upper.append(float(np.log(upper)))
        else:
            opt_p0.append(float(value))
            opt_lower.append(float(lower))
            opt_upper.append(float(upper))

    def unpack(values: np.ndarray) -> list[float]:
        out: list[float] = []
        for name, value in zip(names, values):
            out.append(float(np.exp(value)) if name in log_set else float(value))
        return out

    def residual(values: np.ndarray) -> np.ndarray:
        yhat_local = _as_float_array(model(x_arr, *unpack(values)))
        res = yhat_local - y_arr
        if sigma is not None:
            res = res / _as_float_array(sigma)
        return res

    opt = least_squares(
        residual,
        np.asarray(opt_p0, dtype=float),
        bounds=(np.asarray(opt_lower, dtype=float), np.asarray(opt_upper, dtype=float)),
        x_scale="jac",
        max_nfev=int(max_nfev),
    )
    if not opt.success:
        raise RuntimeError(f"least_squares failed: {opt.message}")

    popt = np.asarray(unpack(opt.x), dtype=float)
    yhat = _as_float_array(model(x_arr, *popt))
    return ParametricFit(
        values={name: float(value) for name, value in zip(names, popt)},
        errors=_least_squares_errors(
            param_names=names,
            popt=popt,
            jacobian=np.asarray(opt.jac, dtype=float),
            residual=np.asarray(opt.fun, dtype=float),
            log_params=log_set,
        ),
        yhat=yhat,
        rmse=rmse(y_arr, yhat),
        chi2=chi2(y_arr, yhat),
        r2=r2_score(y_arr, yhat),
        method=str(method),
        raw_result=opt,
    )


def least_squares_with_standard_errors(
    fun: Callable[[np.ndarray], np.ndarray],
    p0: Sequence[float],
    bounds: tuple[Sequence[float], Sequence[float]],
) -> tuple[np.ndarray, np.ndarray, Any]:
    res = least_squares(fun, x0=p0, bounds=bounds)
    p = np.asarray(res.x, dtype=float)
    se = np.full_like(p, np.nan, dtype=float)
    if res.jac is not None:
        jac = np.asarray(res.jac, dtype=float)
        dof = max(0, jac.shape[0] - jac.shape[1])
        if dof > 0:
            s2 = float(res.cost * 2 / dof)
            try:
                cov = np.linalg.inv(jac.T @ jac) * s2
                se = np.sqrt(np.diag(cov))
            except np.linalg.LinAlgError:
                pass
    return p, se, res


def stderr_from_single_param_jacobian(
    y: np.ndarray,
    yhat: np.ndarray,
    jac: np.ndarray,
    *,
    n_params: int = 1,
) -> float | None:
    y_arr = _as_float_array(y)
    yhat_arr = _as_float_array(yhat)
    jac_arr = _as_float_array(jac)
    if y_arr.size == 0 or yhat_arr.size != y_arr.size or jac_arr.size != y_arr.size:
        return None
    if not np.all(np.isfinite(y_arr)) or not np.all(np.isfinite(yhat_arr)) or not np.all(np.isfinite(jac_arr)):
        return None
    dof = int(y_arr.size) - int(n_params)
    if dof <= 0:
        return None
    sse = float(np.sum((y_arr - yhat_arr) ** 2))
    jtj = float(np.sum(jac_arr ** 2))
    if not np.isfinite(sse) or not np.isfinite(jtj) or jtj <= 0.0:
        return None
    var = sse / (float(dof) * jtj)
    if not np.isfinite(var) or var < 0.0:
        return None
    return float(np.sqrt(var))


def format_value_error(value: object, error: object | None = None, *, unit: str = "") -> str:
    value_f = float(value)
    unit_suffix = f" {unit}" if unit else ""
    if error is None:
        return f"{value_f:.6g}{unit_suffix}"
    error_f = float(error)
    if not np.isfinite(error_f):
        return f"{value_f:.6g}{unit_suffix}"
    return f"{value_f:.6g} +/- {error_f:.3g}{unit_suffix}"
