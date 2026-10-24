.PHONY: demo live setup setup-web test test-engine test-web typecheck build dev-api dev-web pipeline status

PY := .venv/bin/python
CONFIG ?= configs/ar_lbd.yaml
DB ?= data/ar_lbd/web/varfold.sqlite

setup:            ## create the venv and install engine + web dependencies
	python3.12 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e "engine[ml,serve,dev]" -r web/backend/requirements.txt
	cd web/frontend && npm ci

setup-web:        ## lightweight: just what `make demo` needs (no PyTorch, no model)
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -r web/backend/requirements.txt
	cd web/frontend && npm ci

test: test-engine test-web typecheck

test-engine:
	cd engine && ../$(PY) -m pytest -q

test-web:
	cd web/backend && ../../$(PY) -m pytest -q

typecheck:
	cd web/frontend && npx tsc --noEmit

build:            ## production frontend bundle
	cd web/frontend && npm run build

dev-api:          ## API on :8000 against the local bundle
	cd web/backend && VARFOLD_DB=$(abspath $(DB)) ../../.venv/bin/uvicorn app.main:app_factory --factory --reload --port 8000

dev-web:          ## Vite dev server on :5173 (proxies /api to :8000)
	cd web/frontend && npx vite

demo: build      ## site only, from the bundled demo dataset (no model needed) → http://127.0.0.1:8000
	cd web/backend && VARFOLD_DB=$(abspath demo/varfold-demo.sqlite) VARFOLD_STATIC_DIR=$(abspath web/frontend/dist) \
	  ../../.venv/bin/uvicorn app.main:app_factory --factory --port 8000

live:             ## site + live ESMFold folding on this computer → http://127.0.0.1:8000
	./scripts/live-local.sh

pipeline:         ## full research run: curate -> fold -> llr -> validate -> analyze -> export
	.venv/bin/varfold --config $(CONFIG) all

status:
	.venv/bin/varfold --config $(CONFIG) status
