import pytest

from synthetic import build_synthetic_run


@pytest.fixture
def synthetic_run(tmp_path):
    return build_synthetic_run(tmp_path)
