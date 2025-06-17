import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def test_data_path():
    return Path(__file__).parent / "test_data"

@pytest.fixture(scope="session")
def test_data_system_1_path(test_data_path):
    return test_data_path / "system_1"

@pytest.fixture
def test_fchk_path(tmp_path):
    fchk_path = tmp_path / "test.fchk"
    fchk_path.touch()
    return fchk_path

@pytest.fixture
def test_tab_path(tmp_path):
    tab_path = tmp_path / "test.tab"
    tab_path.touch()
    return tab_path