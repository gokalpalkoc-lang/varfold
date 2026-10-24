"""Showcase web app: public read-only API over the research dataset + live-fold proxy."""

from __future__ import annotations

import asyncio
import gzip
import logging
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .db import LABELS, SORTABLE, Dataset
from .live import (LiveError, RateLimiter, WorkerClient, clean_sequence, client_key, point_mutations,
                   superpose_onto)
from .settings import Settings

log = logging.getLogger("varfold.web")


class LiveFoldRequest(BaseModel):
    sequence: str = Field(min_length=1, max_length=5000)


def create_app(settings: Settings | None = None, worker_transport=None) -> FastAPI:
    settings = settings or Settings()

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        app.state.dataset = Dataset(settings.db_path)
        app.state.wt_pdb = gzip.decompress(app.state.dataset.structure_gz("WT")).decode()
        app.state.worker = (WorkerClient(settings.worker_url, settings.worker_token, settings.worker_timeout_s,
                                         transport=worker_transport) if settings.live_enabled else None)
        app.state.limiter = RateLimiter(settings.live_per_hour)
        app.state.live_slots = asyncio.Semaphore(settings.live_concurrency)
        log.info("serving %s (%s); live folding %s", settings.db_path,
                 app.state.dataset.meta["target"]["title"], "on" if app.state.worker else "off")
        yield
        if app.state.worker:
            await app.state.worker.close()

    app = FastAPI(title="VarFold API", version="1.0", lifespan=lifespan, docs_url="/api/docs",
                  openapi_url="/api/openapi.json", redoc_url=None,
                  description="Precomputed ESMFold structures and structural-confidence metrics for "
                              "clinically classified missense variants, plus a rate-limited live-fold endpoint.")

    @app.middleware("http")
    async def headers(request: Request, call_next):
        response = await call_next(request)
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("Referrer-Policy", "strict-origin-when-cross-origin")
        path = request.url.path
        if request.method == "GET" and path.startswith("/api/") and not path.startswith("/api/live"):
            response.headers.setdefault("Cache-Control", "public, max-age=300")
        elif path.startswith("/assets/") and response.status_code == 200:
            # Vite emits content-hashed filenames, so these never change in place.
            response.headers.setdefault("Cache-Control", "public, max-age=31536000, immutable")
        return response

    def ds(request: Request) -> Dataset:
        return request.app.state.dataset

    @app.get("/api/health", include_in_schema=False)
    async def health(request: Request):
        return {"status": "ok", "exported_at": ds(request).meta["exported_at"]}

    @app.get("/api/meta", summary="Target, model, provenance, annotations and summary statistics")
    async def meta(request: Request):
        d = ds(request)
        return {**d.meta, "wild_type": {"sequence": d.wt_sequence, "plddt": d.wt_plddt},
                "live": {"enabled": settings.live_enabled, "max_length": settings.live_max_length}}

    @app.get("/api/variants", summary="List variants with filters, sorting and paging")
    async def variants(request: Request,
                       label: list[str] | None = Query(None, description="repeat or comma-separate"),
                       q: str | None = Query(None, max_length=20, description="variant id prefix or position"),
                       breaking: bool = False,
                       sort: str = "position", order: str = Query("asc", pattern="^(asc|desc)$"),
                       limit: int = Query(50, ge=1, le=500), offset: int = Query(0, ge=0)):
        labels = [x for part in (label or []) for x in part.split(",") if x]
        if bad := set(labels) - LABELS:
            raise HTTPException(422, f"unknown label(s): {sorted(bad)}")
        if sort not in SORTABLE:
            raise HTTPException(422, f"sort must be one of {sorted(SORTABLE)}")
        total, items = ds(request).list_variants(labels, q, breaking, sort, order == "desc", limit, offset)
        return {"total": total, "limit": limit, "offset": offset, "items": items}

    @app.get("/api/variants/{variant_id}", summary="One variant with per-residue pLDDT")
    async def variant(variant_id: str, request: Request):
        v = ds(request).variant(variant_id)
        if v is None:
            raise HTTPException(404, "unknown variant")
        return v

    @app.get("/api/variants/{variant_id}/structure.pdb", summary="Predicted structure (PDB, superposed on WT)",
             response_class=Response)
    async def structure(variant_id: str, request: Request):
        d = ds(request)
        etag = f'"{d.meta["exported_at"]}-{variant_id}"'
        if request.headers.get("if-none-match") == etag:
            return Response(status_code=304, headers={"ETag": etag})
        blob = d.structure_gz(variant_id)
        if blob is None:
            raise HTTPException(404, "unknown variant")
        headers = {"ETag": etag, "Cache-Control": "public, max-age=86400", "Vary": "Accept-Encoding",
                   "Content-Disposition": f'inline; filename="{variant_id}.pdb"'}
        if "gzip" in request.headers.get("accept-encoding", ""):
            return Response(blob, media_type="chemical/x-pdb", headers={**headers, "Content-Encoding": "gzip"})
        return Response(gzip.decompress(blob), media_type="chemical/x-pdb", headers=headers)

    @app.get("/api/positions", summary="Per-residue aggregates and breaking points")
    async def positions(request: Request):
        return ds(request).positions()

    @app.get("/api/live/status", summary="Is the live inference engine reachable?")
    async def live_status(request: Request):
        worker: WorkerClient | None = request.app.state.worker
        if worker is None:
            return {"available": False, "configured": False}
        return {**await worker.status(), "configured": True}

    @app.post("/api/live/fold", summary="Fold a sequence (dataset hit or live GPU inference)")
    async def live_fold(body: LiveFoldRequest, request: Request):
        d = ds(request)
        state = request.app.state
        try:
            seq = clean_sequence(body.sequence, settings.live_max_length)
            mutations = point_mutations(d.wt_sequence, seq, d.meta["construct"]["start"])
            if hit := d.find_sequence(seq):
                return {"source": "dataset", "variant_id": hit, "sequence": seq, "mutations": mutations}
            if state.worker is None:
                raise LiveError(503, "live folding is not enabled on this deployment")
            key = client_key(request.client.host if request.client else None)
            state.limiter.check(key)
            try:
                async with state.live_slots:
                    result = await state.worker.fold(seq)
            except LiveError as exc:
                if exc.status in (503, 504):  # nothing was folded; don't charge the quota
                    state.limiter.refund(key)
                raise
        except LiveError as exc:
            headers = {"Retry-After": str(exc.retry_after)} if exc.retry_after else None
            return JSONResponse({"detail": exc.detail}, status_code=exc.status, headers=headers)
        pdb = superpose_onto(result["pdb"], state.wt_pdb) if mutations is not None else result["pdb"]
        return {"source": "live", "sequence": seq, "mutations": mutations, "pdb": pdb,
                "plddt": result["plddt"], "mean_plddt": result["mean_plddt"], "ptm": result["ptm"],
                "runtime_s": result["runtime_s"], "cached": result.get("cached", False),
                "device": result.get("device"),
                # Live predictions are only directly comparable to the dataset under identical settings.
                "comparable": result.get("fingerprint") == d.meta["model"]["fingerprint"]}

    if settings.static_dir:
        static = settings.static_dir.resolve()
        index = static / "index.html"
        if (static / "assets").is_dir():
            app.mount("/assets", StaticFiles(directory=static / "assets"), name="assets")

        @app.get("/{path:path}", include_in_schema=False)
        async def spa(path: str):
            if path.startswith("api/"):
                raise HTTPException(404)
            candidate = (static / path).resolve()
            if path and candidate.is_file() and candidate.is_relative_to(static):
                return FileResponse(candidate)
            return FileResponse(index, headers={"Cache-Control": "no-cache"})

    return app


def app_factory() -> FastAPI:
    """Entry point for ``uvicorn app.main:app_factory --factory`` (settings from env)."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    return create_app()
