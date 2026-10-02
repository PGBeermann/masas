# Predictor EI-MS — espectros de masas por impacto electrónico a partir de SMILES

Aplicación **Node.js** que recibe una estructura química (dibujada en **JSME** o escrita como SMILES), predice su espectro de masas **EI a 70 eV** y explica, mediante **mecanismos de reacción** (formalismo de sitio de carga/radical de McLafferty), la formación de los iones más importantes.

| Componente | Tecnología |
|---|---|
| Editor molecular | JSME 2024.04.29 (Bienfait & Ertl, *J. Cheminform.* 2013, 5, 24), servido localmente |
| Química (parseo SMILES, aromaticidad, anillos, canonización, dibujo SVG) | RDKit.js (MinimalLib/WASM) 2026.03 |
| Motor de fragmentación | Reglas mecanísticas propias (`src/ms/`) |
| Servidor | Express 5 + Helmet + compresión + caché LRU + límite de tasa |
| Salidas | Espectro (gráfico/tabla), mecanismos ilustrados, CSV, JSON, **JCAMP-DX 5.01**, impresión/PDF |

---

## 1. Instalación local

```bash
npm install
npm start            # http://127.0.0.1:3100
npm test             # validación con 20 compuestos de referencia + 30 de robustez
```

Requiere Node ≥ 18 (probado con Node 22).

## 2. Despliegue en el VPS (Hostinger, Ubuntu)

```bash
# 1) Node LTS y PM2
curl -fsSL https://deb.nodesource.com/setup_22.x | sudo -E bash -
sudo apt-get install -y nodejs nginx
sudo npm i -g pm2

# 2) Código
sudo mkdir -p /var/www/eims && sudo chown $USER /var/www/eims
cd /var/www/eims && unzip ~/eims-predictor.zip && npm ci --omit=dev

# 3) Proceso
pm2 start ecosystem.config.js && pm2 save && pm2 startup   # ejecutar la línea que imprime pm2 startup

# 4) Nginx + TLS
sudo cp deploy/nginx-eims.conf /etc/nginx/sites-available/eims
sudo ln -s /etc/nginx/sites-available/eims /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl reload nginx
sudo certbot --nginx -d espectros.su-dominio.edu.pa
```

Variables de entorno (`ecosystem.config.js`): `PORT` (3100), `HOST` (127.0.0.1, sólo accesible vía Nginx), `BASE_PATH` (para publicar en un subdirectorio, p. ej. `/espectros`), `RATE_LIMIT` (solicitudes/min por IP y por instancia).

## 2-bis. Despliegue con Dokploy (recomendado si el VPS ya usa Dokploy)

Dokploy construye la imagen desde el `Dockerfile` incluido y publica el servicio a través de su Traefik (TLS automático con Let's Encrypt). **No se usan** `ecosystem.config.js` (PM2) ni `deploy/nginx-eims.conf`.

1. Suba el proyecto a un repositorio Git (GitHub/GitLab/Gitea) **sin** `node_modules` (ya está en `.gitignore`).
2. En Dokploy: *Project → Create Service → Application*.
3. *General → Provider*: el repositorio y la rama. *Build Type*: **Dockerfile** (ruta `Dockerfile`, contexto `.`).
4. *Environment* (opcional; el Dockerfile ya trae estos valores por defecto):
   ```
   PORT=3100
   HOST=0.0.0.0
   TRUST_PROXY=1
   RATE_LIMIT=60
   BASE_PATH=
   ```
5. *Domains → Add Domain*: host `espectros.su-dominio.edu.pa`, path `/`, **Container Port 3100**, HTTPS activado, certificado *Let's Encrypt*. El registro DNS A del subdominio debe apuntar a la IP del VPS.
6. *Deploy*. El contenedor queda sano cuando `GET /api/health` responde `{"ok":true}` (HEALTHCHECK del Dockerfile).
7. Escalado: *Advanced → Replicas* (cada réplica es independiente; el límite de tasa se aplica por réplica).

Diferencias respecto de la instalación con PM2/Nginx:

| Aspecto | PM2 + Nginx | Dokploy |
|---|---|---|
| Interfaz de escucha | `HOST=127.0.0.1` | `HOST=0.0.0.0` (obligatorio dentro del contenedor) |
| Proxy de confianza | `TRUST_PROXY=loopback` | `TRUST_PROXY=1` (Traefik es un salto; necesario para que el límite de tasa vea la IP real del cliente) |
| TLS / dominio | Certbot + Nginx | Traefik de Dokploy |
| Reinicio y apagado | PM2 | Docker (`SIGTERM` → cierre ordenado del servidor) |

## 3. API

`POST /api/predict` — cuerpo `{"smiles": "CCCCC(C)=O", "threshold": 0.5}`

Respuesta (extracto):

```json
{
  "smiles": "CCCCC(C)=O", "formula": "C6H12O", "nominalMass": 100, "monoisotopicMass": 100.08882, "rdb": 1,
  "basePeak": { "mz": 43, "formula": "C2H3O" },
  "spectrum": [ { "mz": 58, "intensity": 63.1, "assignments": [ { "formula": "C3H6O", "isotope": false, "share": 1 } ] } ],
  "explanations": [ { "mz": 58, "formula": "C3H6O", "ionType": "OE⁺• (impar-electrónico)", "pathway": [ { "ruleName": "Reordenamiento de McLafferty (cetona)", "mechanism": ["…"], "ionSmiles": "…", "ionSvg": "<svg…>" } ] } ],
  "references": [ … ]
}
```

`GET /api/health` — estado del motor.

## 4. Arquitectura del motor

```
SMILES ─► RDKit (validación, sanitización, aromaticidad, SSSR) ─► MolGraph (átomos pesados + H implícitos)
       ─► reglas primarias sobre M⁺• ─► reglas secundarias (≤ 2 generaciones, regla del electrón par)
       ─► agregación por composición elemental ─► patrón isotópico (convolución) ─► normalización (pico base = 100)
       ─► explicación: ruta mecanística, estructura del ion y del neutro (RDKit SVG), átomos retenidos/enlaces rotos
```

| Archivo | Contenido |
|---|---|
| `src/chem/elements.js` | Masas monoisotópicas (AME 2020) y abundancias (IUPAC 2013); fórmula de Hill, RDB |
| `src/chem/isotopes.js` | Patrón isotópico por convolución polinomial (Cl, Br, S, Si…) |
| `src/chem/molgraph.js` | Grafo molecular desde el JSON de RDKit |
| `src/ms/species.js` | Ion/neutro como subgrafo + transferencias de H, órdenes de enlace, cargas y radicales → SMILES/SVG |
| `src/ms/stability.js` | Escalas de estabilidad de cationes y radicales (regla de Stevenson, pérdida del radical mayor) |
| `src/ms/rules.js` | Reglas primarias |
| `src/ms/secondary.js` | Fragmentaciones consecutivas |
| `src/ms/engine.js` | Orquestación, estabilidad de M⁺•, espectro y explicaciones |

### Reglas implementadas

| # | Mecanismo | Ejemplo diagnóstico |
|---|---|---|
| 1 | α-Escisión iniciada por el sitio radical en N, O, S, X (iones iminio/oxonio/tionio) y pérdida de H• | Trietilamina m/z 86; 2-butanol m/z 45 |
| 2 | α-Escisión del carbonilo (ion acilio) y escisión inductiva (R⁺) | 2-hexanona m/z 43; C₆H₅CO⁺ m/z 105 |
| 3 | Escisión inductiva C–X (haluros, éteres, nitro) | 1-bromobutano m/z 57; nitrobenceno m/z 77 |
| 4 | Ion halonio cíclico de 5 miembros | 1-bromohexano m/z 135/137 |
| 5 | Escisión bencílica → tropilio; escisión alílica | Tolueno/etilbenceno m/z 91 |
| 6 | Escisión σ C–C con preferencia por ramificación | n-Alcanos m/z 43, 57, 71, 85 |
| 7 | Reordenamiento de McLafferty (C=O, C=N, C≡N, C=C, aromático) con ambas retenciones de carga | 2-hexanona m/z 58; ácido butanoico m/z 60 |
| 8 | Retro-Diels–Alder (ciclohexenos, incl. fusionados a arenos) | Ciclohexeno m/z 54 |
| 9 | Apertura de carbociclos (−•CH₃, −C₂H₄) | Ciclohexeno m/z 67; ciclohexano m/z 56 |
| 10 | Pérdidas neutras: H₂O, H₂S, HX, CO y HCO• (fenoles), HCN (anilinas/azinas), NO• (nitroarenos), ceteno (acetanilidas/acetatos de arilo), CH₂O y •CH₃ (anisoles) | Acetanilida m/z 93; fenol m/z 66 |
| S1 | Acilio → R⁺ + CO | 105 → 77 |
| S2 | Reacción del onio (pérdida de alqueno desde iminio/oxonio) | 86 → 58; 59 → 31 |
| S3 | Carbocationes alquilo: −H₂, −CH₄, −alqueno | 43 → 41; 57 → 41 |
| S4 | Degradación de iones aromáticos (−C₂H₂, −CO) | 91 → 65 → 39; 77 → 51; 93 → 65 |

Cada ion se clasifica como **OE⁺•** o **EE⁺** (RDB entero/semientero), se verifica la **regla del nitrógeno** y se informan m/z nominal y exacto (restando la masa del electrón).

## 5. Validación

`npm test` compara la predicción con los picos principales de 20 espectros EI de referencia (NIST WebBook, SRD 69):

* **Pico base correcto: 20/20.**
* **Picos de referencia predichos (≥ 2 %): 80/91 (88 %).**

| Compuesto | Predicho (m/z:%) | Referencia NIST (aprox.) |
|---|---|---|
| 2-Hexanona | 43:100 58:63 100:16 85:9 | 43:100 58:60 71:8 100:8 85:5 |
| Tolueno | 91:100 92:73 65:12 39:7 | 91:100 92:75 65:10 39:8 |
| Butirofenona | 105:100 120:63 77:63 148:41 51:37 | 105:100 120:45 77:40 51:15 148:15 |
| Trietilamina | 86:100 58:62 101:23 30:22 | 86:100 58:30 30:20 101:18 |
| 1-Bromohexano | 135:100 137:97 85:44 43:15 | 135:100 137:95 43:80 41:50 55:40 85:20 |
| Ciclohexeno | 67:100 54:78 82:70 | 67:100 54:60 82:40 41:30 |
| n-Octano | 43:100 57:83 41:63 71:39 85:29 | 43:100 57:35 41:30 85:25 71:20 114:6 |
| Dietil éter | 31:100 59:91 74:66 29:49 | 31:100 59:60 74:60 29:50 45:30 |
| Nitrobenceno | 77:100 123:98 51:58 93:32 | 77:100 123:60 51:55 65:15 93:15 |
| Acetanilida | 93:100 135:33 43:32 | 93:100 135:50 43:30 65:15 66:15 |

### Alcance y limitaciones (declaración metodológica)

1. Las intensidades son **semicuantitativas**: factores heurísticos calibrados contra espectros de referencia, no cálculos cinéticos RRKM/QET. Sirven para docencia, hipótesis de asignación y cribado; **no sustituyen** la comparación con bibliotecas NIST/Wiley ni la medición experimental.
2. No se modelan aún: efecto *orto*, reordenamientos de esqueleto profundos (esteroides, terpenos policíclicos), fragmentación de policiclos fusionados saturados (p. ej., decalina), retro-reacciones de heterociclos (p. ej., pérdida de CH₃NCO en cafeína), derivados TMS específicos, ni iones de carga doble.
3. Las estructuras de iones de reordenamiento se dibujan con su esqueleto neutro y se rotulan ⁺•; los iones descritos sólo por fórmula (p. ej., C₅H₅⁺) no se dibujan.
4. Los átomos se citan con el índice del dibujo RDKit (base 0).

### Cómo añadir una regla

Cree en `src/ms/rules.js` una función `(g) => eventos[]` que devuelva objetos `{ rule, ruleName, ion: Species, neutral, score, steps, refs, oe }` y agréguela a `PRIMARY_RULES`. Las fragmentaciones consecutivas se añaden en `src/ms/secondary.js`. Ejecute `npm test` para verificar que la calibración no se degrada.

## 6. Referencias

1. McLafferty, F. W.; Tureček, F. *Interpretation of Mass Spectra*, 4.ª ed.; University Science Books, 1993.
2. Gross, J. H. *Mass Spectrometry: A Textbook*, 3.ª ed.; Springer, 2017. doi:10.1007/978-3-319-54398-7
3. Silverstein, R. M.; Webster, F. X.; Kiemle, D. J.; Bryce, D. L. *Spectrometric Identification of Organic Compounds*, 8.ª ed.; Wiley, 2014.
4. de Hoffmann, E.; Stroobant, V. *Mass Spectrometry: Principles and Applications*, 3.ª ed.; Wiley, 2007.
5. McLafferty, F. W. Mass spectrometric analysis. Molecular rearrangements. *Anal. Chem.* 1959, 31, 82–87.
6. Stevenson, D. P. Ionization and dissociation by electronic impact. *Discuss. Faraday Soc.* 1951, 10, 35–45.
7. Gohlke, R. S.; McLafferty, F. W. Mass spectrometric analysis. Aliphatic amines. *Anal. Chem.* 1962, 34, 1281–1287.
8. McLafferty, F. W. Mass spectrometric analysis. Aliphatic halogenated compounds. *Anal. Chem.* 1962, 34, 2–15.
9. Rylander, P. N.; Meyerson, S.; Grubb, H. M. Organic ions in the gas phase. II. The tropylium ion. *J. Am. Chem. Soc.* 1957, 79, 842–846.
10. Tureček, F.; Hanuš, V. Retro-Diels–Alder reaction in mass spectrometry. *Mass Spectrom. Rev.* 1984, 3, 85–152.
11. Friedel, R. A.; Shultz, J. L.; Sharkey, A. G. Mass spectra of alcohols. *Anal. Chem.* 1956, 28, 926–934.
12. Karni, M.; Mandelbaum, A. The 'even-electron rule'. *Org. Mass Spectrom.* 1980, 15, 53–64.
13. Meija, J. et al. Isotopic compositions of the elements 2013 (IUPAC Technical Report). *Pure Appl. Chem.* 2016, 88, 293–306.
14. Wang, M. et al. The AME 2020 atomic mass evaluation (II). *Chin. Phys. C* 2021, 45, 030003.
15. Bienfait, B.; Ertl, P. JSME: a free molecule editor in JavaScript. *J. Cheminform.* 2013, 5, 24.
16. McDonald, R. S.; Wilks, P. A. JCAMP-DX: A standard form for exchange of infrared spectra in computer readable form. *Appl. Spectrosc.* 1988, 42, 151–162; Lampen, P. et al. JCAMP-DX for mass spectrometry. *Appl. Spectrosc.* 1994, 48, 1545–1552.
17. NIST Chemistry WebBook, NIST Standard Reference Database Number 69. https://webbook.nist.gov

Licencias de terceros: RDKit (BSD-3), JSME (BSD-3), Express (MIT), Helmet (MIT).
