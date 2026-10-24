#!/usr/bin/env bash
# Run the whole platform on one computer: ESMFold worker + website with live folding.
#   ./scripts/live-local.sh                 # uses your own export if present, else the demo bundle
# Then open http://127.0.0.1:8000  (Ctrl-C stops both processes)
set -euo pipefail
cd "$(dirname "$0")/.."
DB="${VARFOLD_DB:-data/ar_lbd/web/varfold.sqlite}"
[[ -f "$DB" ]] || DB="demo/varfold-demo.sqlite"
[[ -d web/frontend/dist ]] || (cd web/frontend && npm run build)

export VARFOLD_WORKER_TOKEN="${VARFOLD_WORKER_TOKEN:-$(python3 -c 'import secrets; print(secrets.token_hex(32))')}"
.venv/bin/varfold serve --host 127.0.0.1 --port 8100 &
WORKER=$!
trap 'kill $WORKER 2>/dev/null' EXIT

echo "Loading ESMFold (first start downloads ~8 GB of weights)..."
until curl -sf -H "Authorization: Bearer $VARFOLD_WORKER_TOKEN" http://127.0.0.1:8100/health >/dev/null; do
  kill -0 $WORKER 2>/dev/null || { echo "worker failed to start" >&2; exit 1; }
  sleep 2
done

echo "Site: http://127.0.0.1:8000  (dataset: $DB)"
cd web/backend
VARFOLD_DB="$(realpath "../../$DB")" VARFOLD_STATIC_DIR="$(realpath ../frontend/dist)" \
VARFOLD_WORKER_URL=http://127.0.0.1:8100 VARFOLD_WORKER_TOKEN="$VARFOLD_WORKER_TOKEN" VARFOLD_LIVE_PER_HOUR=1000 \
  ../../.venv/bin/uvicorn app.main:app_factory --factory --host 127.0.0.1 --port 8000
