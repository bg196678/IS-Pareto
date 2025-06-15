import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def test_data_path():
    return Path(__file__).parent / "test_data"

@pytest.fixture(scope="session")
def test_data_system_1_path(test_data_path):
    return test_data_path / "system_1"
