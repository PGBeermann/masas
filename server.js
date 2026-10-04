'use strict';
/**
 * EI-MS Predictor — servidor HTTP (Express).
 *   POST /api/predict  { smiles, threshold? }  → espectro + explicación mecanística
 *   GET  /api/health
 * Diseñado para ejecutarse detrás de Nginx (proxy inverso) con PM2 en un VPS.
 */
const path = require('path');
const express = require('express');
const helmet = require('helmet');
const compression = require('compression');
const initRDKit = require('@rdkit/rdkit');
const { predict } = require('./src/ms/engine');

const PORT = Number(process.env.PORT || 3100);
const HOST = process.env.HOST || '127.0.0.1';
const BASE = (process.env.BASE_PATH || '').replace(/\/$/, ''); // p. ej. "/espectros" si se publica en un subdirectorio

const app = express();
// Detrás de Nginx local: 'loopback'. Detrás de Traefik (Dokploy/Docker): 1 salto → TRUST_PROXY=1
const TP = process.env.TRUST_PROXY || 'loopback';
app.set('trust proxy', /^\d+$/.test(TP) ? Number(TP) : TP === 'true' ? true : TP);
app.disable('x-powered-by');
app.use(helmet({
  contentSecurityPolicy: {
    useDefaults: true,
    directives: {
      // JSME es una aplicación GWT: inyecta scripts y requiere eval / estilos en línea.
      'script-src': ["'self'", "'unsafe-inline'", "'unsafe-eval'"],
      'style-src': ["'self'", "'unsafe-inline'", 'https://fonts.googleapis.com'],
      'font-src': ["'self'", 'https://fonts.gstatic.com', 'data:'],
      'img-src': ["'self'", 'data:', 'blob:'],
      'frame-src': ["'self'"],
      'connect-src': ["'self'"]
    }
  },
  crossOriginEmbedderPolicy: false
}));
app.use(compression());
app.use(express.json({ limit: '16kb' }));

/* --- limitador de tasa simple en memoria (por IP) --- */
const WINDOW_MS = 60_000; const MAX_REQ = Number(process.env.RATE_LIMIT || 60);
const hits = new Map();
function rateLimit(req, res, next) {
  const now = Date.now(); const ip = req.ip;
  const rec = hits.get(ip) || { n: 0, t: now };
  if (now - rec.t > WINDOW_MS) { rec.n = 0; rec.t = now; }
  rec.n++; hits.set(ip, rec);
  if (rec.n > MAX_REQ) return res.status(429).json({ error: 'Demasiadas solicitudes. Intente en un minuto.' });
  next();
}
setInterval(() => { const now = Date.now(); for (const [k, v] of hits) if (now - v.t > WINDOW_MS) hits.delete(k); }, WINDOW_MS).unref();

/* --- caché LRU de resultados --- */
const CACHE_MAX = 300; const cache = new Map();
function cacheGet(k) { if (!cache.has(k)) return null; const v = cache.get(k); cache.delete(k); cache.set(k, v); return v; }
function cacheSet(k, v) { cache.set(k, v); if (cache.size > CACHE_MAX) cache.delete(cache.keys().next().value); }

let RDKit = null;
const router = express.Router();

router.get('/api/health', (req, res) => res.json({ ok: !!RDKit, rdkit: RDKit ? RDKit.version() : null, uptime: process.uptime() }));

router.post('/api/predict', rateLimit, (req, res) => {
  if (!RDKit) return res.status(503).json({ error: 'El motor químico aún se está inicializando.' });
  const smiles = String((req.body && req.body.smiles) || '').trim();
  const threshold = Math.min(10, Math.max(0, Number(req.body && req.body.threshold) || 0.5));
  if (!smiles) return res.status(400).json({ error: 'Debe proporcionar un SMILES.' });
  if (smiles.length > 500 || /[^A-Za-z0-9@+\-\[\]\(\)=#$%\/\\.:*~]/.test(smiles)) return res.status(400).json({ error: 'SMILES con caracteres no permitidos o demasiado largo.' });
  const key = smiles + '|' + threshold;
  const hit = cacheGet(key); if (hit) return res.json(hit);
  try {
    const out = predict(RDKit, smiles, { threshold });
    cacheSet(key, out);
    res.json(out);
  } catch (e) {
    res.status(422).json({ error: e.message || 'No fue posible procesar la estructura.' });
  }
});

router.use('/jsme', express.static(path.join(__dirname, 'node_modules', 'jsme-editor'), { maxAge: '30d' }));
// index.html sin caché (siempre se revalida); app.js y style.css se versionan con ?v= en index.html,
// de modo que cada actualización del programa llega al navegador sin necesidad de forzar la recarga.
router.use(express.static(path.join(__dirname, 'public'), {
  maxAge: '1h',
  setHeaders: (res, file) => { if (file.endsWith('.html')) res.setHeader('Cache-Control', 'no-cache'); }
}));

app.use(BASE || '/', router);
app.use((err, req, res, next) => { // eslint-disable-line no-unused-vars
  if (err.type === 'entity.too.large') return res.status(413).json({ error: 'Solicitud demasiado grande.' });
  console.error(err); res.status(500).json({ error: 'Error interno.' });
});

initRDKit().then(R => {
  RDKit = R;
  const server = app.listen(PORT, HOST, () => console.log(`EI-MS Predictor en http://${HOST}:${PORT}${BASE || '/'}  (RDKit ${R.version()})`));
  // Apagado ordenado (docker stop / redeploy de Dokploy envía SIGTERM)
  for (const sig of ['SIGTERM', 'SIGINT']) process.on(sig, () => { server.close(() => process.exit(0)); setTimeout(() => process.exit(0), 5000).unref(); });
}).catch(e => { console.error('No se pudo inicializar RDKit:', e); process.exit(1); });
