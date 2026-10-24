import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// In development the API runs separately (uvicorn on :8000); Vite proxies /api to it.
export default defineConfig({
  plugins: [react()],
  server: {
    proxy: { "/api": { target: process.env.VARFOLD_API ?? "http://127.0.0.1:8000", changeOrigin: true } },
  },
  build: {
    target: "es2022",
    chunkSizeWarningLimit: 1200, // 3Dmol is large; it is split out and loaded lazily
  },
});
