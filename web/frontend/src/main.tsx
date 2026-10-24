import { StrictMode, type ReactNode } from "react";
import { createRoot } from "react-dom/client";
import { BrowserRouter, Link, Route, Routes } from "react-router";
import { api, useAsync } from "./lib/api";
import { ErrorNote, Layout, Loading } from "./components/ui";
import { Overview } from "./pages/Overview";
import { Variants } from "./pages/Variants";
import { VariantDetail } from "./pages/VariantDetail";
import { LiveFold } from "./pages/LiveFold";
import { Methods } from "./pages/Methods";
import "./styles.css";

function App() {
  const meta = useAsync(api.meta, "meta");
  const page = (render: (m: NonNullable<typeof meta.data>) => ReactNode) =>
    meta.data ? render(meta.data) : meta.error ? <ErrorNote error={meta.error} /> : <Loading what="dataset" />;
  return (
    <Routes>
      <Route element={<Layout meta={meta.data} />}>
        <Route index element={page((m) => <Overview meta={m} />)} />
        <Route path="variants" element={page((m) => <Variants meta={m} />)} />
        <Route path="variants/:id" element={page((m) => <VariantDetail meta={m} />)} />
        <Route path="live" element={page((m) => <LiveFold meta={m} />)} />
        <Route path="methods" element={page((m) => <Methods meta={m} />)} />
        <Route path="*" element={<div className="stack-sm"><h1>Page not found</h1><Link to="/">Back to the overview</Link></div>} />
      </Route>
    </Routes>
  );
}

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <BrowserRouter>
      <App />
    </BrowserRouter>
  </StrictMode>,
);
