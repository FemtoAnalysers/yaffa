import os
import subprocess
from pathlib import Path
import pytest
from dotenv import load_dotenv


ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "build"
BINARY = BUILD_DIR / "bin" / "SimulateSource"
ROOTDIFF = ROOT / "scripts" / "rootdiff.py"


def run(cmd, cwd=None):
    """Run command and fail test on error."""
    print(f"\n[RUN] {' '.join(map(str, cmd))}")
    subprocess.run(cmd, cwd=cwd, check=True)


@pytest.fixture(scope="session", autouse=True)
def setup_environment():
    # Load .env (like `source ../.env`)
    env_path = ROOT / ".env"
    if env_path.exists():
        load_dotenv(env_path)

    # Ensure build dir exists
    BUILD_DIR.mkdir(parents=True, exist_ok=True)

    # Configure + build
    run(["cmake", "-DLOG_LEVEL=40", str(ROOT)], cwd=BUILD_DIR)
    run(["make"], cwd=BUILD_DIR)

    # Ensure binary exists
    assert BINARY.exists(), f"Missing binary: {BINARY}"


def run_simulation(cfg, out_root, ref_root):
    run([str(BINARY), cfg])
    run(["python3", str(ROOTDIFF), out_root, ref_root])


# TODO: fix this test
@pytest.mark.skip(reason="Unable to import libCATS.so")
def test_ceca():
    run_simulation(
        "cfg_test_ceca.yaml",
        "source.root",
        "source_ref.root",
    )


# TODO: fix this test
@pytest.mark.skip(reason="Unable to import libCATS.so")
def test_ceca_3b():
    run_simulation(
        "cfg_test_ceca_3b.yaml",
        "source_3b.root",
        "source_3b_ref.root",
    )


# TODO: fix this test
@pytest.mark.skip(reason="Unable to import libCATS.so")
def test_cecapaper():
    run_simulation(
        "cfg_test_cecapaper.yaml",
        "source_cecapaper.root",
        "source_cecapaper_ref.root",
    )


# TODO: fix this test
@pytest.mark.skip(reason="Unable to import libCATS.so")
def test_cecapaper_3b():
    run_simulation(
        "cfg_test_cecapaper_3b.yaml",
        "source_cecapaper_3b.root",
        "source_cecapaper_3b_ref.root",
    )
