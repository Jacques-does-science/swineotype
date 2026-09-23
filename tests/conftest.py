"""Fixtures shared across the swineotype test suite.

Constants and builders live in tests/helpers.py, imported as a top-level
module: pytest prepends tests/ to sys.path, and `tests` is a name already
taken by an unrelated package in some environments.
"""
import sys
from pathlib import Path

import pytest

# scripts/resolver_refs.py is the one FASTA reader/writer and manifest checker;
# the golden tests and the BLAST fixtures both import it.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))


@pytest.fixture
def patched_stages(monkeypatch):
    """Patch process_one's collaborators. Yields a setter."""
    from swineotype import main as main_mod

    def setup(s1, s2=None, forbid_stage2=False):
        monkeypatch.setattr(main_mod, "stage1_score", lambda *a, **k: s1)
        monkeypatch.setattr(main_mod, "ensure_unix_line_endings", lambda p, t: p)
        if forbid_stage2:
            monkeypatch.setattr(main_mod, "stage2_resolver_call",
                                lambda *a, **k: pytest.fail("Stage 2 must not run here"))
        else:
            monkeypatch.setattr(main_mod, "stage2_resolver_call", lambda *a, **k: s2)
        return main_mod

    return setup
