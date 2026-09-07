#!/usr/bin/env python3
import json
import platform
import sys
from pathlib import Path

import scipy
import sklearn

root = Path(__file__).resolve().parents[1]
payload = {
    "python": sys.version,
    "python_executable": sys.executable,
    "platform": platform.platform(),
    "scikit_learn": sklearn.__version__,
    "scipy": scipy.__version__,
    "phase1_benchmark_runtime_dependencies": "Python standard library only",
    "validation_dependency": "scikit-learn metrics",
}
path = root / "results" / "logs" / "phase1_python_environment.json"
path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
