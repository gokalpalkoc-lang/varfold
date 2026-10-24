export function num(v: number | null | undefined, digits = 1, signed = false): string {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  const s = v.toFixed(digits);
  return signed && v > 0 ? `+${s}` : s.replace("-", "−");
}

export function pval(v: number | null | undefined): string {
  if (v === null || v === undefined) return "—";
  if (v < 1e-3) {
    const [m, e] = v.toExponential(1).split("e");
    return `${m}×10${superscript(e)}`;
  }
  return v.toFixed(3);
}

function superscript(exp: string): string {
  const map: Record<string, string> = { "-": "⁻", "+": "", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
  return [...exp].map((c) => map[c] ?? c).join("");
}

export function pct(v: number | null | undefined, digits = 0): string {
  return v === null || v === undefined ? "—" : `${(v * 100).toFixed(digits)}%`;
}

export function date(iso: string | undefined): string {
  if (!iso) return "—";
  return new Date(iso).toLocaleDateString(undefined, { year: "numeric", month: "short", day: "numeric" });
}
