from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run_script(idx: int, total: int, script: Path):
    print(f"[{idx}/{total}] Running {script.name} ...")
    subprocess.run([sys.executable, str(script)], check=True)


def main():
    project_root = Path(__file__).resolve().parent
    scripts = [
        project_root / "scripts" / "dynamic_gravity.py",
        project_root / "scripts" / "geodesic_light.py",
        project_root / "scripts" / "cosmic_expansion.py",
        project_root / "scripts" / "run_quantum_clock.py",
        project_root / "scripts" / "thermodynamics.py",
        project_root / "scripts" / "holography.py",
        project_root / "scripts" / "unification.py",
    ]

    total = len(scripts)
    for i, script in enumerate(scripts, start=1):
        run_script(i, total, script)

    print("=== Orchestration complete ===")
    # Show verdict from unification
    final_report = project_root / "outputs" / "THEORY_OF_EVERYTHING.md"
    if final_report.exists():
        print(f"Final report: {final_report}")
        print(final_report.read_text(encoding="utf-8").splitlines()[0:7])


if __name__ == "__main__":
    main()
