"""Input-path configuration for the pipeline.

All locations of the (large, public) input files are read from a JSON file, by default
`config_paths.json` next to this script, or from the file named by the environment variable
GYN_FIDELITY_PATHS. Copy `config_paths.example.json`, edit the paths to where the public
files were downloaded, and save it as `config_paths.json`. Relative paths are resolved
against the directory of the JSON file. No path is hard-coded in the analysis scripts.
"""
import json, os
from pathlib import Path

PIPE = Path(__file__).resolve().parent
_cfg_file = Path(os.environ.get('GYN_FIDELITY_PATHS', PIPE / 'config_paths.json'))
if not _cfg_file.exists():
    raise SystemExit(f'Path configuration not found: {_cfg_file}\n'
                     f'Copy config_paths.example.json to config_paths.json and edit it.')
_cfg = json.load(open(_cfg_file, encoding='utf-8'))


def path(key):
    """Return the configured path for `key` (absolute, or relative to the config file)."""
    if key not in _cfg:
        raise KeyError(f'"{key}" missing from {_cfg_file}')
    p = Path(_cfg[key]).expanduser()
    return p if p.is_absolute() else (_cfg_file.parent / p).resolve()


def has(key):
    return key in _cfg and bool(_cfg[key])


def paths(key):
    """List-valued entry (e.g. tcga_file_dirs)."""
    v = _cfg.get(key, [])
    v = v if isinstance(v, list) else [v]
    out = []
    for x in v:
        p = Path(x).expanduser()
        out.append(p if p.is_absolute() else (_cfg_file.parent / p).resolve())
    return out


def locate_tcga(local_path):
    """Resolve a GDC STAR file named <file_id>__<file_name> (manifest column local_path holds the file name)
    by searching the configured tcga_file_dirs; returns the first existing path, or the path in the first dir."""
    name = Path(str(local_path)).name
    dirs = paths('tcga_file_dirs')
    for d in dirs:
        if (d / name).exists():
            return d / name
    return dirs[0] / name if dirs else Path(name)
