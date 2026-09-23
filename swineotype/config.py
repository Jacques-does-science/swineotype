import os
import tempfile
from importlib.resources import files
from pathlib import Path

import yaml

# --- Default Configuration ---

DEFAULT_CONFIG = {
    "wzxwzy_fasta": "suis_wzxwzy_whitelist.fasta",
    "resolver_refs_fasta": "suis_resolver_refs.fasta",
    "tmp_dir": "",  # empty = derive per-run (see load_config)
    "plurality": 0.60,
    "delta": 100,
    # Only call a type whose serotype-specific gene (wzy) was found.
    # wzx is conserved across serotypes and cannot carry a call on its own.
    "require_wzy": 1,
    "min_pid": 85.0,
    "min_cov": 0.80,
    "min_res_pid": 90.0,
    "min_res_alen": 300,
    "keep_debug": 1,
    "gzip_debug": 0,
    # Species this tool assigns serotypes for. A cps reference tagged with any
    # other species identifies a different organism, not an S. suis serotype.
    "target_species": "Streptococcus suis",
    # The resolvable families. A type in neither is its own family and needs
    # no within-family resolution.
    "pair_1_14": {"1", "14"},
    "pair_2_1_2": {"2", "1/2"},
}

# --- Configuration Loading ---

def data_dir() -> Path:
    """The reference data: shipped inside the package, or $SWINEOTYPE_HOME/data."""
    if "SWINEOTYPE_HOME" in os.environ:
        return Path(os.environ["SWINEOTYPE_HOME"]) / "data"
    return Path(str(files("swineotype") / "data"))


def load_config(config_file: str | None = None) -> dict:
    """
    Loads configuration from a YAML file, filling in with defaults.
    """
    config = DEFAULT_CONFIG.copy()

    if config_file:
        with open(config_file, "r") as f:
            user_config = yaml.safe_load(f)
            if user_config:
                config.update(user_config)

    # --- Environment Variable Overrides ---

    # Coerce against the DEFAULT's type. `type(value)(raw)` was wrong for the
    # set-valued keys -- set("1,14") yields {'1', ',', '4'} -- and cannot
    # express a "not set" sentinel at all.
    for key, default in DEFAULT_CONFIG.items():
        env_var = f"SWINEO_{key.upper()}"
        if env_var not in os.environ:
            continue
        raw = os.environ[env_var]
        if isinstance(default, (set, frozenset)):
            config[key] = {t.strip() for t in raw.split(",") if t.strip()}
        elif isinstance(default, int):
            config[key] = int(raw)
        elif isinstance(default, float):
            config[key] = float(raw)
        else:
            config[key] = raw

    # --- Path Resolution ---

    config["data_dir"] = data_dir()
    config["wzxwzy_fasta"] = config["data_dir"] / config["wzxwzy_fasta"]
    config["resolver_refs_fasta"] = config["data_dir"] / config["resolver_refs_fasta"]

    # tmp_dir: an explicit setting (config file or SWINEO_TMP_DIR) wins and is
    # flagged so the CLI does not override it. Otherwise fall back to the
    # system temp dir -- NOT the install tree, which may be read-only and is
    # shared between unrelated runs. The CLI normally replaces this with
    # <out_dir>/.swineotype_cache.
    explicit = bool(config.get("tmp_dir"))
    config["tmp_dir_explicit"] = explicit
    config["tmp_dir"] = Path(config["tmp_dir"]) if explicit \
        else Path(tempfile.gettempdir()) / "swineotype_cache"
    config["tmp_dir"].mkdir(parents=True, exist_ok=True)

    return config
