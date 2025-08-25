from swineotype.config import load_config

def test_load_config_defaults():
    config = load_config()
    assert config["plurality"] == 0.60
    assert config["delta"] == 100
    assert config["min_pid"] == 85.0
