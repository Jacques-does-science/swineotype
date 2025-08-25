from swineotype.config import load_config

from unittest.mock import patch


def test_load_config_defaults():
    config = load_config()
    assert config["plurality"] == 0.60
    assert config["delta"] == 100
    assert config["min_pid"] == 85.0

def test_load_config_env_vars():
    with patch.dict("os.environ", {"SWINEO_PLURALITY": "0.75", "SWINEO_DELTA": "200"}):
        config = load_config()
        assert config["plurality"] == 0.75
        assert config["delta"] == 200

