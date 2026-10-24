"""Runtime configuration, all from environment variables (see deploy/.env.example)."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


def _int(name: str, default: int) -> int:
    return int(os.environ.get(name, default))


@dataclass
class Settings:
    db_path: Path = field(default_factory=lambda: Path(os.environ.get("VARFOLD_DB", "/data/varfold.sqlite")))
    static_dir: Path | None = field(default_factory=lambda: Path(p) if (p := os.environ.get("VARFOLD_STATIC_DIR")) else None)
    # Live inference worker on the research PC, reached over Tailscale.
    worker_url: str | None = field(default_factory=lambda: os.environ.get("VARFOLD_WORKER_URL") or None)
    worker_token: str = field(default_factory=lambda: os.environ.get("VARFOLD_WORKER_TOKEN", ""))
    worker_timeout_s: int = field(default_factory=lambda: _int("VARFOLD_WORKER_TIMEOUT", 300))
    live_max_length: int = field(default_factory=lambda: _int("VARFOLD_LIVE_MAX_LENGTH", 400))
    # Per-client budget for live folds that miss the precomputed dataset.
    live_per_hour: int = field(default_factory=lambda: _int("VARFOLD_LIVE_PER_HOUR", 6))
    live_concurrency: int = field(default_factory=lambda: _int("VARFOLD_LIVE_CONCURRENCY", 4))

    @property
    def live_enabled(self) -> bool:
        return bool(self.worker_url and self.worker_token)
