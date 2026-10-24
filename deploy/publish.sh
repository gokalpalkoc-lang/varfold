#!/usr/bin/env bash
# Upload a freshly exported dataset bundle to the server and restart the app.
#   VARFOLD_SERVER=deploy@your-server ./deploy/publish.sh [path/to/varfold.sqlite]
set -euo pipefail
: "${VARFOLD_SERVER:?set VARFOLD_SERVER=user@host}"
REMOTE_DIR="${VARFOLD_REMOTE_DIR:-/opt/varfold/deploy}"
DB="${1:-data/ar_lbd/web/varfold.sqlite}"
[[ -f "$DB" ]] || { echo "no bundle at $DB — run: varfold export" >&2; exit 1; }

rsync -az --progress "$DB" "$VARFOLD_SERVER:$REMOTE_DIR/data/varfold.sqlite.new"
# Atomic swap; the app opens the file read-only/immutable, so restart it.
ssh "$VARFOLD_SERVER" "cd '$REMOTE_DIR' && mv data/varfold.sqlite.new data/varfold.sqlite && docker compose restart app"
echo "published $(basename "$DB") to $VARFOLD_SERVER"
