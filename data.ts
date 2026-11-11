// Synthetic FIRMS-like reference data, progressive sampling scenarios and the two measures.
export const K = 10, B = 12, MINN = 5, N = 5000;
export const FR = [.02, .04, .07, .11, .16, .22, .3, .4, .52, .65];
// Square grid in Web Mercator metres.
export const CELL = 400000, X0 = -13950000, YT = 6400000, GX = 17, GY = 9, NC = GX * GY;
const RE = 6378137, D2R = Math.PI / 180;
export const merc = (lng: number, lat: number) => ({ x: RE * lng * D2R, y: RE * Math.log(Math.tan(Math.PI / 4 + lat * D2R / 2)) });
export const unmerc = (x: number, y: number) => ({ lng: x / RE / D2R, lat: (2 * Math.atan(Math.exp(y / RE)) - Math.PI / 2) / D2R });
export const cellAt = (lat: number, lng: number) => {
  const m = merc(lng, lat), c = Math.floor((m.x - X0) / CELL), r = Math.floor((YT - m.y) / CELL);
  return c < 0 || r < 0 || c >= GX || r >= GY ? -1 : r * GX + c;
};
const OUT: [number, number][] = [[-124.7,48.4],[-124.1,46.2],[-124.5,42.8],[-124.4,40.4],[-122.5,37.8],[-121.9,36.6],[-120.6,34.5],[-118.5,34],[-117.1,32.5],[-114.7,32.7],[-111,31.3],[-108.2,31.3],[-106.5,31.8],[-104.5,29.6],[-103,29],[-101.4,29.8],[-99.5,27.5],[-97.4,25.9],[-97.2,27.8],[-94.7,29.3],[-93.8,29.7],[-91,29.2],[-89.4,29],[-89.2,30.2],[-87.5,30.3],[-85.3,29.7],[-84,30],[-82.8,28.9],[-82.6,27.5],[-81.8,26.1],[-81,25.1],[-80.4,25.3],[-80.1,26.8],[-80.6,28.4],[-81.4,30.6],[-81,32],[-79,33.6],[-77.9,34],[-76.5,34.7],[-75.5,35.3],[-76,36.9],[-75.1,38.4],[-74,39.5],[-74,40.5],[-72,41],[-70,41.7],[-70.7,42.7],[-70.2,43.7],[-68.5,44.3],[-67,44.8],[-67.8,47.1],[-69.2,47.4],[-70.8,45.3],[-74.8,45],[-79,43.3],[-83,42],[-82.4,45.3],[-84.7,46.5],[-89.5,48],[-95.2,49],[-123,49],[-122.8,48.3],[-124.7,48.4]];
function inUS(lo: number, la: number) {
  let c = false;
  for (let i = 0, j = OUT.length - 1; i < OUT.length; j = i++) {
    const a = OUT[i], b = OUT[j];
    if ((a[1] > la) !== (b[1] > la) && lo < (b[0] - a[0]) * (la - a[1]) / (b[1] - a[1]) + a[0]) c = !c;
  }
  return c;
}
function rng(s: number) { return () => { s |= 0; s = s + 0x6D2B79F5 | 0; let t = Math.imul(s ^ s >>> 15, 1 | s); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
const R = rng(42), gs = () => Math.sqrt(-2 * Math.log(R() + 1e-12)) * Math.cos(6.2832 * R());
// lng, lat, spread (deg), mean ln(FRP), share
const CL = [[-121.8,40.3,1.3,3.3,.16],[-117.3,34.1,1,3,.1],[-114.5,45.5,1.8,2.7,.12],[-106.5,39.5,1.6,2.6,.08],[-110.5,34.2,1.6,2.9,.1],[-98,31.3,1.8,2.5,.12],[-96.5,37.5,1.5,2.3,.1],[-83.5,32.5,1.8,2.4,.1],[-81.5,28.5,1,2.2,.06],[0,0,0,2.2,.08]];
export interface Pt { lat: number; lng: number; lf: number; c: number; r: number; cell: number; b: number }
export const P: Pt[] = [];
while (P.length < N) {
  let r = R(), c = 0;
  while (c < 9 && r > CL[c][4]) { r -= CL[c][4]; c++; }
  const q = CL[c]; let lo: number, la: number;
  if (c === 9) { lo = -125 + R() * 58; la = 25 + R() * 24; } else { lo = q[0] + q[2] * gs(); la = q[1] + q[2] * .8 * gs(); }
  if (!inUS(lo, la)) continue;
  const m = merc(lo, la), cx = Math.floor((m.x - X0) / CELL), cy = Math.floor((YT - m.y) / CELL);
  if (cx < 0 || cy < 0 || cx >= GX || cy >= GY) continue;
  const lf = q[3] + .65 * gs();
  P.push({ lat: la, lng: lo, lf, c: cx, r: cy, cell: cy * GX + cx, b: Math.max(0, Math.min(B - 1, Math.floor((lf - .5) / 5.5 * B))) });
}
export const refN = new Int32Array(NC);
export const refH = Array.from({ length: NC }, () => new Float64Array(B));
P.forEach(p => { refN[p.cell]++; refH[p.cell][p.b]++; });
refH.forEach((h, i) => { if (refN[i]) for (let b = 0; b < B; b++) h[b] /= refN[i]; });
const Q = rng(7), U = P.map(() => Q());
function makeOrder(sc: string) {
  const s = sc === "structural" || sc === "combined", d = sc === "distributional" || sc === "combined";
  const wc = (p: Pt) => s ? .25 + 2.75 * (1 - p.c / (GX - 1)) : 1;             // west favoured, east starved
  const wf = (p: Pt) => d && p.r < GY / 2 ? Math.exp(-1.2 * Math.max(-2, Math.min(3, p.lf - 2.5))) : 1; // north favours low FRP
  const key = P.map((p, i) => -Math.log(U[i]) / wf(p) / wc(p)), idx = P.map((_, i) => i);
  if (s) return idx.sort((a, b) => key[a] - key[b]);
  const gr: number[][] = Array.from({ length: NC }, () => []); idx.forEach(i => gr[P[i].cell].push(i));
  const f = new Float64Array(N);
  gr.forEach(a => { a.sort((i, j) => key[i] - key[j]); a.forEach((i, k) => f[i] = (k + .5) / a.length + 1e-6 * U[i]); });
  return idx.sort((a, b) => f[a] - f[b]);
}
function jsd(p: ArrayLike<number>, q: ArrayLike<number>) {
  let s = 0;
  for (let i = 0; i < B; i++) { const a = p[i], b = q[i], m = (a + b) / 2; if (a > 0) s += a * Math.log2(a / m) / 2; if (b > 0) s += b * Math.log2(b / m) / 2; }
  return s;
}
export interface Stage { ns: number; n: Int32Array; hh: Int32Array; st: (number | null)[]; js: (number | null)[] }
const cache: Record<string, { ord: number[]; M: Stage[] }> = {};
export function data(sc: string) {
  if (cache[sc]) return cache[sc];
  const ord = makeOrder(sc), cnt = new Int32Array(NC), h = new Int32Array(NC * B), M: Stage[] = []; let pos = 0;
  for (let k = 0; k < K; k++) {
    const ns = Math.round(N * FR[k]);
    while (pos < ns) { const p = P[ord[pos++]]; cnt[p.cell]++; h[p.cell * B + p.b]++; }
    const n = cnt.slice(), hh = h.slice(), st: (number | null)[] = [], js: (number | null)[] = [];
    for (let c = 0; c < NC; c++) {
      if (!refN[c]) { st.push(null); js.push(null); continue; }
      st.push(n[c] ? Math.max(-1, Math.min(1, Math.log2((n[c] / ns) / (refN[c] / N)) / 2)) : -1);   // structural score
      if (n[c] >= MINN) { const q: number[] = []; for (let b = 0; b < B; b++) q.push(hh[c * B + b] / n[c]); js.push(jsd(refH[c], q)); } else js.push(null); // JSD
    }
    M.push({ ns, n, hh, st, js });
  }
  return cache[sc] = { ord, M };
}
/** 12 log-FRP bins grouped into low / medium / high classes (shares). */
export const frpClasses = (a: ArrayLike<number>) => [[0, 4], [4, 7], [7, 12]].map(([s, e]) => { let t = 0; for (let i = s; i < e; i++) t += a[i]; return t; });
