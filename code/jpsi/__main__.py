"""三个位置参数的命令行入口；结果目录唯一，NPZ 不使用 pickle。"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile
import numpy as np
from .solver import solve


def _flatten(value, prefix=""):
    # 双下划线键名在不使用 pickle 的 NPZ 中保留嵌套结果结构。
    out = {}
    for key, item in value.items():
        name = f"{prefix}__{key}" if prefix else key
        if isinstance(item, dict):
            out.update(_flatten(item, name))
        elif item is not None:
            out[name] = np.asarray(item)
    return out


def _json_finite(value):
    if isinstance(value, dict):
        return {k: _json_finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_finite(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_result(result, output_root):
    """创建唯一结果子目录，绝不静默覆盖已有文件。"""
    # 写入唯一目录，避免新归档覆盖已有数值结果。
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = Path(tempfile.mkdtemp(prefix=f"jpsi_{stamp}_", dir=root))
    with (dest / "result.npz").open("xb") as stream:
        np.savez_compressed(stream, **_flatten(result))
    summary = dict(result_directory=str(dest), metadata=result["metadata"], n_borel=len(result["M2"]), channels={})
    for channel, r in result["channels"].items():
        d = r["alpha_diagnostics"]
        summary["channels"][channel] = dict(alpha=r["alpha_star"], residues=r["residues"].tolist(),
            f_GeV=r["f"].tolist(), weighted_residual=r["weighted_residual"],
            relative_weighted_residual=r["relative_weighted_residual"], source=d["source"],
            fallback=d["fallback"], on_boundary=d["on_boundary"],
            small_alpha_peak_covered=d["small_alpha_peak_covered"])
    summary = _json_finite(summary)
    with (dest / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return summary


def main():
    # 命令行参数依次为 M2_min、M2_max 和 Lambda，单位均为 GeV^2。
    parser = argparse.ArgumentParser(description="J/psi n=0 spectral inversion (GeV^2 inputs)")
    parser.add_argument("M2_min", type=float)
    parser.add_argument("M2_max", type=float)
    parser.add_argument("Lambda", type=float)
    parser.add_argument("--channel", choices=("1delta", "2delta"), default=None)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    try:
        result = solve(args.M2_min, args.M2_max, args.Lambda, channel=args.channel)
        summary = save_result(result, args.output)
    except (ValueError, RuntimeError, FloatingPointError, OSError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(summary, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
