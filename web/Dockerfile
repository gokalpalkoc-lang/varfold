# Showcase app: built React frontend served by the FastAPI backend.
# Build context: web/   (docker compose -f deploy/docker-compose.yml build)

FROM node:24-bookworm-slim AS frontend
WORKDIR /src
COPY frontend/package.json frontend/package-lock.json ./
RUN npm ci --no-audit --no-fund
COPY frontend/ ./
RUN npm run build

FROM python:3.12-slim AS app
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    VARFOLD_STATIC_DIR=/app/static \
    VARFOLD_DB=/data/varfold.sqlite
WORKDIR /app
COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY backend/app ./app
COPY --from=frontend /src/dist ./static
RUN useradd --system --uid 10001 --no-create-home varfold
USER varfold
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/api/health', timeout=3)"
# Only Caddy can reach this port (it is not published), so trusting its
# X-Forwarded-For is safe; the live-fold rate limit keys on the real client IP.
CMD ["uvicorn", "app.main:app_factory", "--factory", "--host", "0.0.0.0", "--port", "8000", \
     "--proxy-headers", "--forwarded-allow-ips", "*"]
