import { useEffect, useMemo, useRef, useState, type KeyboardEvent, type PointerEvent } from "react";

export interface LineSeries {
  key: string;
  label: string;
  values: (number | null)[];
  color: string;
  width?: number;
}

export interface StemSeries {
  label: string;
  values: (number | null)[];
  color: (v: number, index: number) => string;
}

export interface RugRow {
  label: string;
  positions: number[];
  color: string;
}

export interface TooltipRow {
  label: string;
  value: string;
  color?: string;
}

interface Props {
  title: string;
  /** Residue number of index 0 (UniProt numbering). */
  start: number;
  sequence?: string;
  lines?: LineSeries[];
  stems?: StemSeries;
  yDomain: [number, number];
  yTicks?: number[];
  yLabel: string;
  height?: number;
  bands?: { from: number; to: number; color: string }[];
  marker?: number;
  windowRange?: [number, number];
  rug?: RugRow[];
  hover?: number | null;
  onHover?: (position: number | null) => void;
  onSelect?: (position: number) => void;
  extraRows?: (position: number) => TooltipRow[];
  format?: (v: number) => string;
}

const M = { top: 10, right: 12, bottom: 30, left: 46 };
const RUG_ROW = 11;

/** Smallest "nice" step giving at most ``target`` ticks, but never fewer than 3 when possible. */
function niceStep(span: number, target: number): number {
  const steps = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000];
  const max = Math.max(3, target);
  return steps.find((s) => span / s <= max) ?? steps[steps.length - 1];
}

/** Per-residue chart: lines and/or stems over residue position, with a snapping crosshair. */
export function ResidueChart({
  title, start, sequence, lines = [], stems, yDomain, yTicks, yLabel, height = 220, bands, marker, windowRange,
  rug = [], hover: hoverProp, onHover, onSelect, extraRows, format = (v) => v.toFixed(1),
}: Props) {
  const wrap = useRef<HTMLDivElement>(null);
  const [width, setWidth] = useState(720);
  const [localHover, setLocalHover] = useState<number | null>(null);
  const [pointerX, setPointerX] = useState<number | null>(null);
  const hover = hoverProp !== undefined ? hoverProp : localHover;
  const setHover = (p: number | null) => (onHover ? onHover(p) : setLocalHover(p));

  useEffect(() => {
    if (!wrap.current) return;
    const ro = new ResizeObserver(([entry]) => setWidth(Math.max(280, entry.contentRect.width)));
    ro.observe(wrap.current);
    return () => ro.disconnect();
  }, []);

  const n = lines[0]?.values.length ?? stems?.values.length ?? 0;
  const rugH = rug.length * RUG_ROW + (rug.length ? 6 : 0);
  const plotW = width - M.left - M.right;
  const plotTop = M.top + rugH;
  const plotH = height - plotTop - M.bottom;
  const band = plotW / Math.max(1, n);
  const x = (i: number) => M.left + (i + 0.5) * band;
  const [y0, y1] = yDomain;
  const y = (v: number) => plotTop + plotH - ((Math.max(y0, Math.min(y1, v)) - y0) / (y1 - y0)) * plotH;

  const xTicks = useMemo(() => {
    const step = niceStep(n, Math.floor(plotW / 60));
    const first = Math.ceil(start / step) * step;
    const out: number[] = [];
    for (let p = first; p < start + n; p += step) out.push(p);
    return out;
  }, [n, plotW, start]);
  const ticks = useMemo(() => {
    if (yTicks) return yTicks;
    const step = niceStep(y1 - y0, 4);
    const out: number[] = [];
    for (let v = Math.ceil(y0 / step) * step; v <= y1 + 1e-9; v += step) out.push(+v.toFixed(6));
    return out;
  }, [yTicks, y0, y1]);

  const paths = lines.map((s) => {
    let d = "";
    let pen = false;
    s.values.forEach((v, i) => {
      if (v === null || v === undefined) {
        pen = false;
        return;
      }
      d += `${pen ? "L" : "M"}${x(i).toFixed(1)},${y(v).toFixed(1)}`;
      pen = true;
    });
    return { ...s, d };
  });

  const toPosition = (clientX: number) => {
    const rect = wrap.current!.getBoundingClientRect();
    const px = clientX - rect.left;
    const i = Math.floor((px - M.left) / band);
    return { i: Math.max(0, Math.min(n - 1, i)), px };
  };
  const onMove = (e: PointerEvent) => {
    const { i, px } = toPosition(e.clientX);
    setPointerX(px);
    setHover(start + i);
  };
  const onKey = (e: KeyboardEvent) => {
    const cur = hover ?? start;
    if (e.key === "ArrowRight" || e.key === "ArrowLeft") {
      e.preventDefault();
      const next = Math.max(start, Math.min(start + n - 1, cur + (e.key === "ArrowRight" ? 1 : -1)));
      setPointerX(x(next - start));
      setHover(next);
    } else if (e.key === "Enter" && hover !== null && onSelect) {
      onSelect(hover);
    } else if (e.key === "Escape") {
      setHover(null);
    }
  };

  const hi = hover !== null && hover >= start && hover < start + n ? hover - start : null;
  const stemW = Math.max(1, Math.min(4, band - 1));
  const rows: TooltipRow[] = hi === null ? [] : [
    ...lines.map((s) => ({ label: s.label, value: s.values[hi] == null ? "—" : format(s.values[hi]!), color: s.color })),
    ...(stems ? [{ label: stems.label, value: stems.values[hi] == null ? "—" : format(stems.values[hi]!) }] : []),
    ...(extraRows ? extraRows(start + hi) : []),
  ];
  const tipX = pointerX ?? (hi !== null ? x(hi) : 0);
  const tipLeft = tipX > width - 200 ? tipX - 190 : tipX + 14;

  return (
    <div className="chart" ref={wrap}>
      <svg width={width} height={height} role="img" aria-label={title}>
        {bands?.map((b) => (
          <rect key={b.from} x={M.left} width={plotW} y={y(b.to)} height={y(b.from) - y(b.to)} fill={b.color} opacity={0.09} />
        ))}
        {windowRange && (() => {
          // Clip to the construct, like the engine's window mean.
          const a = Math.max(start, windowRange[0]);
          const b = Math.min(start + n - 1, windowRange[1]);
          return b >= a && (
            <rect x={x(a - start) - band / 2} width={(b - a + 1) * band} y={plotTop} height={plotH} fill="var(--ink)" opacity={0.05} />
          );
        })()}
        <g className="grid">
          {ticks.map((t) => <line key={t} x1={M.left} x2={M.left + plotW} y1={y(t)} y2={y(t)} />)}
        </g>
        {ticks.map((t) => (
          <text key={t} x={M.left - 8} y={y(t)} dy="0.32em" textAnchor="end">{format(t).replace(/\.0$/, "")}</text>
        ))}
        <text className="axis-label" transform={`translate(12 ${plotTop + plotH / 2}) rotate(-90)`} textAnchor="middle">{yLabel}</text>
        {y0 < 0 && y1 > 0 && <line className="baseline" x1={M.left} x2={M.left + plotW} y1={y(0)} y2={y(0)} />}
        <line className="baseline" x1={M.left} x2={M.left + plotW} y1={plotTop + plotH} y2={plotTop + plotH} />
        {xTicks.map((p) => (
          <text key={p} x={x(p - start)} y={plotTop + plotH + 18} textAnchor="middle">{p}</text>
        ))}

        {rug.map((r, k) => (
          <g key={r.label}>
            {r.positions.filter((p) => p >= start && p < start + n).map((p) => (
              <rect key={p} x={x(p - start) - Math.max(1, stemW / 2)} y={M.top + k * RUG_ROW} width={Math.max(2, stemW)}
                    height={RUG_ROW - 3} rx={1} fill={r.color} />
            ))}
          </g>
        ))}

        {stems?.values.map((v, i) => v === null || v === undefined || v === 0 ? null : (
          <rect key={i} x={x(i) - stemW / 2} width={stemW} y={Math.min(y(v), y(0))} height={Math.max(1, Math.abs(y(v) - y(0)))}
                rx={Math.min(2, stemW / 2)} fill={stems.color(v, i)} />
        ))}
        {paths.map((s) => (
          <path key={s.key} d={s.d} fill="none" stroke={s.color} strokeWidth={s.width ?? 2} strokeLinejoin="round" strokeLinecap="round" />
        ))}
        {marker !== undefined && marker >= start && marker < start + n && (
          <g>
            <line x1={x(marker - start)} x2={x(marker - start)} y1={plotTop} y2={plotTop + plotH} stroke="var(--pathogenic)" strokeWidth={1.5} />
            <path d={`M${x(marker - start) - 5},${plotTop - 1} h10 l-5,7 z`} fill="var(--pathogenic)" />
          </g>
        )}
        {hi !== null && (
          <g pointerEvents="none">
            <line x1={x(hi)} x2={x(hi)} y1={plotTop} y2={plotTop + plotH} stroke="var(--ink-2)" strokeWidth={1} />
            {paths.map((s) => s.values[hi] == null ? null : (
              <circle key={s.key} cx={x(hi)} cy={y(s.values[hi]!)} r={4} fill={s.color} stroke="var(--surface)" strokeWidth={2} />
            ))}
          </g>
        )}
        <rect
          x={M.left} y={M.top} width={plotW} height={plotTop + plotH - M.top} fill="transparent" tabIndex={0}
          style={{ cursor: onSelect ? "pointer" : "crosshair", outline: "none" }}
          aria-label={`${title}. Use arrow keys to inspect residues${onSelect ? ", Enter to open" : ""}.`}
          onPointerMove={onMove} onPointerLeave={() => { setHover(null); setPointerX(null); }}
          onClick={(e) => onSelect && onSelect(start + toPosition(e.clientX).i)}
          onKeyDown={onKey} onBlur={() => setHover(null)}
        />
      </svg>
      {hi !== null && (
        <div className="tooltip" style={{ left: tipLeft, top: plotTop + 4 }}>
          <div className="t-head">{sequence?.[hi] ?? ""}{start + hi}</div>
          {rows.map((r) => (
            <div className="t-row" key={r.label}>
              <span className="k">{r.color && <span className="t-key" style={{ background: r.color }} />}{r.label}</span>
              <span className="v">{r.value}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
