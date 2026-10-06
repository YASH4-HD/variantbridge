"""Smoke tests: the Streamlit UI renders without exceptions in each mode.

These check that the app runs, not that any science is right (science is
tested in the other modules). Demo modes are SIMULATED.
"""
from pathlib import Path

import pytest

st_testing = pytest.importorskip("streamlit.testing.v1")
APP = str(Path(__file__).resolve().parents[1] / "app.py")


def _run(branch, source):
    at = st_testing.AppTest.from_file(APP, default_timeout=180).run()
    at.radio[0].set_value(branch).run()
    at.radio[1].set_value(source).run()
    return at


def test_cohort_demo_renders():
    at = _run("CohortOmics - Human cohorts", "Demo data (SIMULATED)")
    assert not at.exception, [e.value for e in at.exception]


def test_amr_demo_renders():
    at = _run("EvoResist-AI - Pathogen/AMR", "Demo data (SIMULATED)")
    assert not at.exception, [e.value for e in at.exception]


def test_cohort_validated_artifacts_without_artifacts_does_not_crash():
    at = _run("CohortOmics - Human cohorts", "Validated artifacts (real data)")
    assert not at.exception, [e.value for e in at.exception]
