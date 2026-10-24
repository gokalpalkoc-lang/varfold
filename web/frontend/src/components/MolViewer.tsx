import { useEffect, useRef, useState } from "react";
import type { GLViewer } from "3dmol";
import { cssVar, deltaColor, plddtColor } from "../lib/colors";

type ThreeDmol = typeof import("3dmol");
let loader: Promise<ThreeDmol> | null = null;
/** 3Dmol is ~1 MB, so it is fetched only when a viewer is first shown. */
function load3Dmol(): Promise<ThreeDmol> {
  loader ??= import("3dmol");
  return loader;
}

export interface MolViewerProps {
  pdb: string | undefined;
  /** Per-residue values indexed by PDB residue number − 1; falls back to the B-factor (pLDDT). */
  values?: number[];
  colorBy?: "plddt" | "delta";
  /** Optional reference model drawn as a translucent ghost (already superposed). */
  overlayPdb?: string;
  /** PDB residue number to show as sticks with a label. */
  highlightResi?: number;
  highlightLabel?: string;
  /** Offset to convert PDB residue numbers into displayed (UniProt) numbering. */
  numberOffset?: number;
  sequence?: string;
  height?: number;
  label: string;
}

export function MolViewer({
  pdb, values, colorBy = "plddt", overlayPdb, highlightResi, highlightLabel, numberOffset = 0, sequence,
  height = 420, label,
}: MolViewerProps) {
  const host = useRef<HTMLDivElement>(null);
  const viewer = useRef<GLViewer | null>(null);
  const [ready, setReady] = useState(false);
  const [failed, setFailed] = useState<string | null>(null);
  const [hover, setHover] = useState<string | null>(null);
  const lastPdb = useRef<string | undefined>(undefined);

  // Create the WebGL viewer once.
  useEffect(() => {
    let disposed = false;
    let observer: ResizeObserver | undefined;
    let themeObserver: MutationObserver | undefined;
    let media: MediaQueryList | undefined;
    let unpaint: (() => void) | undefined;
    load3Dmol()
      .then(($3Dmol) => {
        if (disposed || !host.current) return;
        const v = $3Dmol.createViewer(host.current, { antialias: true });
        if (!v) throw new Error("WebGL unavailable");
        v.setHoverDuration(40);
        viewer.current = v;
        // A transparent canvas isn't honoured everywhere, so paint the theme surface explicitly
        // and repaint whenever the theme (toggle or OS setting) changes.
        const paint = () => {
          v.setBackgroundColor(cssVar("--surface", "#fcfcfb"), 1);
          v.render();
        };
        paint();
        themeObserver = new MutationObserver(paint);
        themeObserver.observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
        media = window.matchMedia("(prefers-color-scheme: dark)");
        media.addEventListener("change", paint);
        unpaint = () => media?.removeEventListener("change", paint);
        observer = new ResizeObserver(() => v.resize());
        observer.observe(host.current);
        setReady(true);
      })
      .catch((e: Error) => setFailed(e.message || "3D viewer failed to load"));
    return () => {
      disposed = true;
      observer?.disconnect();
      themeObserver?.disconnect();
      unpaint?.();
      viewer.current?.clear();
      viewer.current = null;
    };
  }, []);

  // (Re)draw whenever the inputs change; keep the camera unless the structure itself changed.
  useEffect(() => {
    const v = viewer.current;
    if (!ready || !v || !pdb) return;
    const keepView = lastPdb.current === pdb;
    const view = keepView ? v.getView() : null;
    v.removeAllModels();
    v.removeAllLabels();
    if (overlayPdb) {
      const ghost = v.addModel(overlayPdb, "pdb");
      ghost.setStyle({}, { cartoon: { color: "#9a988f", opacity: 0.32 } });
    }
    const model = v.addModel(pdb, "pdb");
    const colorFor = (atom: { resi: number; b: number }) => {
      const val = values?.[atom.resi - 1];
      if (colorBy === "delta") return deltaColor(val ?? 0);
      return plddtColor(val ?? atom.b);
    };
    model.setStyle({}, { cartoon: { colorfunc: colorFor, thickness: 0.4, arrows: true } });
    if (highlightResi !== undefined) {
      model.setStyle({ resi: highlightResi }, {
        cartoon: { colorfunc: colorFor, thickness: 0.4, arrows: true },
        stick: { colorscheme: "orangeCarbon", radius: 0.28 },
      });
      if (highlightLabel) {
        v.addLabel(highlightLabel, {
          fontSize: 13, fontColor: "#ffffff", backgroundColor: "#0b0b0b", backgroundOpacity: 0.82,
          borderRadius: 6, padding: 3, inFront: true, alignment: "bottomCenter",
        } as never, { resi: highlightResi, atom: "CA" });
      }
    }
    model.setHoverable({}, true,
      (atom: { resi: number; resn: string; b: number }) => {
        const pos = atom.resi + numberOffset;
        const aa = sequence?.[atom.resi - 1] ?? atom.resn;
        const val = values?.[atom.resi - 1] ?? atom.b;
        const metric = colorBy === "delta" ? `ΔpLDDT ${val > 0 ? "+" : ""}${val.toFixed(1)}` : `pLDDT ${val.toFixed(1)}`;
        setHover(`${aa}${pos} · ${metric}`);
      },
      () => setHover(null));
    if (view) v.setView(view);
    else v.zoomTo();
    v.render();
    lastPdb.current = pdb;
  }, [ready, pdb, values, colorBy, overlayPdb, highlightResi, highlightLabel, numberOffset, sequence]);

  return (
    <div className="viewer" style={{ height }} role="img" aria-label={label}>
      <div ref={host} style={{ position: "absolute", inset: 0 }} />
      {(!ready || !pdb) && !failed && <div className="viewer-loading">Loading structure…</div>}
      {failed && <div className="viewer-loading error">{failed}</div>}
      <div className="viewer-overlay">
        <span className="viewer-hint">{hover ?? "Drag to rotate · scroll to zoom"}</span>
      </div>
    </div>
  );
}
