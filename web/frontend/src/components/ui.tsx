import { useEffect, useState, type ReactNode } from "react";
import { NavLink, Outlet, Link } from "react-router";
import type { Label, Meta } from "../lib/api";
import { LABEL_NAMES, labelVar, PLDDT_BANDS } from "../lib/colors";
import { date } from "../lib/format";

export function Layout({ meta }: { meta: Meta | undefined }) {
  return (
    <div className="shell">
      <a className="sr-only" href="#main">Skip to content</a>
      <header className="topbar">
        <div className="container topbar-inner">
          <Link to="/" className="brand" aria-label="VarFold home">
            <img src="/favicon.svg" alt="" />
            <span>VarFold</span>
          </Link>
          <nav className="nav" aria-label="Primary">
            <NavLink to="/" end>Overview</NavLink>
            <NavLink to="/variants">Variants</NavLink>
            <NavLink to="/live">Live fold</NavLink>
            <NavLink to="/methods">Methods</NavLink>
          </nav>
          <ThemeToggle />
        </div>
      </header>
      <main id="main" className="container">
        <Outlet />
      </main>
      <footer className="footer">
        <div className="container">
          <span>
            {meta ? <>{meta.target.gene} · UniProt <a href={`https://www.uniprot.org/uniprotkb/${meta.target.uniprot}`} target="_blank" rel="noreferrer">{meta.target.uniprot}</a> · dataset exported {date(meta.exported_at)}</> : "VarFold"}
          </span>
          <span>
            Predictions: ESMFold (Lin et al., 2023) · <a href="/api/docs">API</a> · Research use only — not for clinical decisions
          </span>
        </div>
      </footer>
    </div>
  );
}

function ThemeToggle() {
  const [theme, setTheme] = useState<string | undefined>(() => document.documentElement.dataset.theme);
  useEffect(() => {
    if (theme) document.documentElement.dataset.theme = theme;
    else delete document.documentElement.dataset.theme;
    try {
      if (theme) localStorage.setItem("varfold-theme", theme);
      else localStorage.removeItem("varfold-theme");
    } catch {
      /* storage unavailable: theme still applies for this visit */
    }
  }, [theme]);
  const systemDark = window.matchMedia?.("(prefers-color-scheme: dark)").matches;
  const isDark = theme ? theme === "dark" : systemDark;
  return (
    <button className="btn icon" onClick={() => setTheme(isDark ? "light" : "dark")}
            aria-label={`Switch to ${isDark ? "light" : "dark"} theme`} title="Toggle theme">
      <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" aria-hidden="true">
        {isDark
          ? <><circle cx="12" cy="12" r="4" /><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4" /></>
          : <path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z" />}
      </svg>
    </button>
  );
}

export function StatTile({ label, value, sub }: { label: string; value: ReactNode; sub?: ReactNode }) {
  return (
    <div className="tile">
      <div className="label">{label}</div>
      <div className="value">{value}</div>
      {sub && <div className="sub">{sub}</div>}
    </div>
  );
}

export function ClassBadge({ label }: { label: Label }) {
  return (
    <span className="badge">
      <span className="dot" style={{ background: labelVar(label) }} />
      {LABEL_NAMES[label]}
    </span>
  );
}

export function BreakingBadge({ q }: { q?: number | null }) {
  return <span className="badge break" title={q != null ? `BH q = ${q.toPrecision(2)}` : undefined}>Breaking point</span>;
}

export function PlddtLegend() {
  return (
    <div className="legend" aria-label="pLDDT colour key">
      {PLDDT_BANDS.map((b) => (
        <span className="key" key={b.label}><span className="swatch" style={{ background: b.color }} />{b.label} ({b.range})</span>
      ))}
    </div>
  );
}

export function DeltaLegend() {
  return (
    <div className="legend" aria-label="ΔpLDDT colour key">
      <span className="key"><span className="swatch" style={{ background: "#d03b3b" }} />Confidence lost (−15 or less)</span>
      <span className="key"><span className="swatch" style={{ background: "#d9d8d3" }} />No change</span>
      <span className="key"><span className="swatch" style={{ background: "#256abf" }} />Confidence gained (+15)</span>
    </div>
  );
}

export function LineLegend({ items }: { items: { label: string; color: string; kind?: "line" | "swatch" }[] }) {
  return (
    <div className="legend">
      {items.map((i) => (
        <span className="key" key={i.label}>
          <span className={i.kind === "swatch" ? "swatch" : "line"} style={{ background: i.color }} />{i.label}
        </span>
      ))}
    </div>
  );
}

export function Loading({ what }: { what: string }) {
  return <div className="skeleton" role="status">Loading {what}…</div>;
}

export function ErrorNote({ error }: { error: Error }) {
  return <div className="callout warn" role="alert">Could not load data: {error.message}</div>;
}
