import threading
import time
from types import SimpleNamespace

import numpy as np
import pytest
from fastapi.testclient import TestClient

from varfold.worker import create_app

TOKEN = "s" * 24
AUTH = {"Authorization": f"Bearer {TOKEN}"}


class FakePredictor:
    fingerprint = "fake|recycles=4|fp32"
    device_name = "fake"

    def __init__(self, delay: float = 0.0):
        self.calls = 0
        self.delay = delay

    def predict(self, seq):
        self.calls += 1
        time.sleep(self.delay)
        return SimpleNamespace(pdb="ATOM\n", plddt=np.full(len(seq), 80.0), mean_plddt=80.0, ptm=0.8, runtime_s=0.01)


def make_client(tmp_path, predictor, **kw):
    return TestClient(create_app(lambda: predictor, token=TOKEN, cache_dir=tmp_path / "cache", **kw))


def test_requires_token(tmp_path):
    with pytest.raises(ValueError):
        create_app(FakePredictor, token="short")
    with make_client(tmp_path, FakePredictor()) as c:
        assert c.get("/health").status_code == 401
        assert c.post("/fold", json={"sequence": "MKT"}, headers={"Authorization": "Bearer wrong"}).status_code == 401
        assert c.get("/health", headers=AUTH).json()["status"] == "ok"


def test_fold_validates_and_caches(tmp_path):
    p = FakePredictor()
    with make_client(tmp_path, p, max_length=10) as c:
        assert c.post("/fold", json={"sequence": "MKTX"}, headers=AUTH).status_code == 422
        assert c.post("/fold", json={"sequence": "M" * 11}, headers=AUTH).status_code == 422
        first = c.post("/fold", json={"sequence": "mkt ay"}, headers=AUTH).json()
        second = c.post("/fold", json={"sequence": "MKTAY"}, headers=AUTH).json()
    assert first["cached"] is False and second["cached"] is True
    assert first["length"] == 5 and p.calls == 1


def test_busy_queue_rejects(tmp_path):
    p = FakePredictor(delay=0.5)
    with make_client(tmp_path, p, max_queue=0) as c:
        results = {}

        def slow():
            results["first"] = c.post("/fold", json={"sequence": "MKTAYIAK"}, headers=AUTH).status_code

        t = threading.Thread(target=slow)
        t.start()
        time.sleep(0.15)
        # The GPU is taken and the queue holds nobody else: reject immediately.
        busy = c.post("/fold", json={"sequence": "MKTAYIAKQ"}, headers=AUTH)
        t.join()
    assert results["first"] == 200
    assert busy.status_code == 503 and busy.headers["retry-after"] == "30"
