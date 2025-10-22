import os
from pathlib import Path
import yaml
import site

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

def get_root_dir() -> Path:
    """
    Finds the project root directory.

    The logic is as follows:
    1. If SWINEOTYPE_HOME is set, use it.
    2. Check for a `data` directory next to the `swineotype` package.
       This works for editable installs (`pip install -e .`).
    3. Check for a `data` directory in `sys.prefix` for standard installs.
    """
    if "SWINEOTYPE_HOME" in os.environ:
        return Path(os.environ["SWINEOTYPE_HOME"])

    package_path = Path(__file__).parent

    # Editable install: <root>/swineotype
    editable_install_data_path = package_path.parent / "data"
    if editable_install_data_path.exists():
        return package_path.parent

    # Standard install: <prefix>/lib/pythonX.Y/site-packages/swineotype
    # and <prefix>/share/swineotype/data
    for sp_path in site.getsitepackages():
        if package_path.is_relative_to(sp_path):
            prefix = Path(sp_path).parent.parent.parent
            share_data_path = prefix / "share" / "swineotype" / "data"
            if share_data_path.exists():
                return prefix / "share" / "swineotype"

    raise FileNotFoundError("Could not locate the `data` directory. "
                            "Please set the SWINEOTYPE_HOME environment variable.")


def load_config(config_file: str | None = None) -> dict:
    """
    Loads configuration from a YAML file, filling in with defaults.
    """
    config = DEFAULT_CONFIG.copy()
    root_dir = get_root_dir()

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

    config["data_dir"] = root_dir / "data"
    config["wzxwzy_fasta"] = config["data_dir"] / config["wzxwzy_fasta"]
    config["resolver_refs_fasta"] = config["data_dir"] / config["resolver_refs_fasta"]
    config["tmp_dir"] = config["data_dir"] / "tmp"
    config["tmp_dir"].mkdir(parents=True, exist_ok=True)

    return config
