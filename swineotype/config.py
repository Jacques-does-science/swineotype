import os
from pathlib import Path
import yaml

# --- Default Configuration ---

DEFAULT_CONFIG = {
    "data_dir": "data",
    "wzxwzy_fasta": "suis_wzxwzy_whitelist.fasta",
    "resolver_refs_fasta": "suis_resolver_refs.fasta",
    "tmp_dir": "results/db_cache",
    "plurality": 0.60,
    "delta": 100,
    "require_agreement": 1,
    "min_pid": 85.0,
    "min_cov": 0.80,
    "min_res_pid": 90.0,
    "min_res_alen": 300,
    "keep_debug": 1,
    "gzip_debug": 0,
    "clean_temp": 0,
    "ambig_set": {"1", "14", "2", "1/2"},
    "pair_1_14": {"1", "14"},
    "pair_2_1_2": {"2", "1/2"},
}

# --- Configuration Loading ---

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

    for key, value in config.items():
        env_var = f"SWINEO_{key.upper()}"
        if env_var in os.environ:
            config[key] = type(value)(os.environ[env_var])

    # --- Path Resolution ---

    config["data_dir"] = Path(config["data_dir"]).resolve()
    config["wzxwzy_fasta"] = config["data_dir"] / config["wzxwzy_fasta"]
    config["resolver_refs_fasta"] = config["data_dir"] / config["resolver_refs_fasta"]
    config["tmp_dir"] = Path(config["tmp_dir"]).resolve()
    config["tmp_dir"].mkdir(parents=True, exist_ok=True)

    return config
