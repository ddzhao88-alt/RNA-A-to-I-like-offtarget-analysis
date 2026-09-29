from pathlib import Path
import yaml


def load_config(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return cfg


def resolve(root, value):
    p = Path(value)
    return p if p.is_absolute() else root / p
