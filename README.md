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
| `src/ms/rules.js` | Reglas primarias (1–11) |
| `src/ms/rules_ext.js` | Reglas ampliadas v1.1 (12–15): anillos/ion distónico, retro-reacción de imidas, efecto orto, pérdida de ROH |
| `src/ms/rules_beauchamp.js` | Reglas v1.2 (B1–B7) tomadas de Beauchamp, *Basics of Mass Spectroscopy* (MS_chapter.pdf): migración de enlaces múltiples, alquinos, alcanos −RH, eliminaciones, fragmentación alquénica consecutiva |
| `src/ms/secondary.js` | Fragmentaciones consecutivas |
| `src/ms/arrows.js` | Dibujo de mecanismos (v1.3): flechas curvas de 1 e⁻ (media cabeza) y 2 e⁻ (completa), carga ⊕ y electrón desapareado • sobre el SVG de RDKit; estado de transición cíclico plegado (McLafferty, eliminaciones, halonio) |
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
| 7 | Reordenamiento de McLafferty (C=O, C=N, C≡N, C=C, aromático) con ambas retenciones de carga; en ésteres, también desde la cadena alcoxílica | 2-hexanona m/z 58; ácido butanoico m/z 60 |
| 7b | **McLafferty + 1** (doble transferencia de H en ésteres → RC(OH)₂⁺) *(v1.1)* | Acetato de etilo m/z 61; benzoato de etilo m/z 123 |
| 7c | **Escisión bencílica/alílica con transferencia de H desde N–H u O–H** (McLafferty aromático con γ-heteroátomo) *(v1.1)* | 2-Feniletanol m/z 92; dopamina m/z 124 |
| 8 | Retro-Diels–Alder (ciclohexenos, incl. fusionados a arenos) | Ciclohexeno m/z 54 |
| 9 | Apertura de carbociclos (−•CH₃, −C₂H₄) | Ciclohexeno m/z 67; ciclohexano m/z 56 |
| 10 | Pérdidas neutras: H₂O, H₂S, HX, CO y HCO• (fenoles), HCN (anilinas/azinas), NO• (nitroarenos), ceteno (acetanilidas/acetatos de arilo), CH₂O y •CH₃ (anisoles) | Acetanilida m/z 93; fenol m/z 66 |
| S1 | Acilio → R⁺ + CO | 105 → 77 |
| S2 | Reacción del onio (pérdida de alqueno desde iminio/oxonio) | 86 → 58; 59 → 31 |
| S3 | Carbocationes alquilo: −H₂, −CH₄, −alqueno | 43 → 41; 57 → 41 |
| S4 | Degradación de iones aromáticos (−C₂H₂, −CO) | 91 → 65 → 39; 77 → 51; 93 → 65 |
| 12 | **α-Escisión en anillos con ion distónico** *(v1.1)*: (A) β-escisión con pérdida de alqueno; (B) transferencia 1,4/1,5-H + β-escisión con pérdida de radical; (C) en éteres cíclicos, expulsión de R₂C=O con carga en el hidrocarburo | Ciclohexanona 55 y 70 → 42; ciclohexanol 57; piperidina 56/57; pirrolidina 43; THF 42 |
| 13 | **Retro-reacción de imidas/ureas cíclicas** (pérdida de R–N=C=O) + −CO y −HCN consecutivos *(v1.1)* | Cafeína 194 → 137 → 109 → 82 → 55 |
| 14 | **Efecto orto** *(v1.1)*: eliminación de ROH/H₂O (donador OH, NH, SH o CH orto a COOR/COOH) y de •OH en o-nitroarenos; se aplica también a iones secundarios | Salicilato de metilo 152 → 120 → 92; aspirina 180 → 138 → 120; o-nitrotolueno 137 → 120 |
| 15 | Pérdida de ROH en ésteres (ion ceteno) *(v1.1)* | [M−ROH]⁺• débil |
| S5 | **Reacción del onio desde el carbono** (R–CH=OH⁺ → CH₂=OH⁺ + alqueno) *(v1.1)* | 2-Butanol 59 → 31 |
| S6 | Cadenas de pérdidas tras reordenamientos: −CO del acilio distónico, −CO tras efecto orto, −•CH₃ del cicloalqueno de deshidratación, −H• de iones CₙH₂ₙ⁺• | Ciclohexanona 70 → 42; ciclohexanol 82 → 67 |

#### Reglas incorporadas en v1.2 (Beauchamp, *Basics of Mass Spectroscopy*, MS_chapter.pdf)

Se revisaron todos los patrones del capítulo frente a la v1.1; sólo se añadieron los ausentes. Los que ya existían (α-escisión en heteroátomos y carbonilos, −CO del acilio, McLafferty en C=O/C=C/C≡N/arenos, RDA, tropilio, halonio de 5 miembros, −H₂O/−HX, α-escisión anular de cicloalcanoles y cicloalquilaminas, pérdida de alqueno/RCOOH en ésteres) no se duplicaron.

| # | Mecanismo | Sección del capítulo | Ejemplo diagnóstico |
|---|---|---|---|
| B1 | Migración del doble/triple enlace en M⁺• (hasta 2 desplazamientos) seguida de escisión alílica/propargílica y McLafferty del isómero | Alquenos §2–4; alquinos §2 | (E)-2-hepteno m/z 69, 70, 41 |
| B2 | Escisión propargílica (catión propargilo ⇄ ciclopropenilo; serie 39, 53, 67), pérdida de H• propargílico [M−1]⁺ y M⁺• débil en alquinos terminales | Alquinos §1, §3 | 1-heptino m/z 39, 95 |
| B3 | McLafferty en alquinos (radical-catión aleno, masa par) | Alquinos §5 | 1-heptino m/z 40 |
| B4 | Alcanos ramificados: eliminación de R–H → ion alqueno CₙH₂ₙ⁺• | Alcanos §8 (p. 19) | 3,4-dimetilhexano m/z 56 |
| B5 | Eliminación de ROH (éteres), RSH (sulfuros), NH₃/RNH₂ (aminas) | Hoja resumen (p. 44) | Dipropil éter m/z 42 |
| B6 | Fragmentación alquénica consecutiva de los iones alqueno formados por −H₂O, −H₂S, −HX, −ROH, −RSH, −RH y McLafferty | Alcoholes §2, tioles §2 | 1-hexanotiol m/z 56, 55 |
| B7 | −(H₂O + C₂H₄) / −(H₂S + C₂H₄) en 1-alcanoles y 1-alcanotioles ≥ C₅ | Alcoholes §2; tioles (ejemplo) | 1-hexanol m/z 56 |
| B8 | Iones onio cíclicos de 5 miembros ampliados a I y S, e iones de 3 miembros (halogenonio/tiiranio) | Halógenos §5; tioles §4 | 1-yodohexano m/z 183; 1-hexanotiol m/z 89, 61 |
| B9 | Doble McLafferty en cetonas con H γ en ambas cadenas | Carbonilos (4-octanona) | 4-octanona m/z 58 |
| B10 | Cationes alilo CₙH₂ₙ₋₁⁺ → −H₂ (41 → 39, 55 → 53, 69 → 67) | Alcanos §7; alquinos §3 | m/z 39 |
| B12 | McLafferty de alquenos con reparto de carga por la regla de Stevenson entre los dos alquenos (EI estimada por tamaño y grado de sustitución del C=C); antes la carga en el alqueno expulsado se trataba como minoritaria | Alquenos §4 (hept-1-eno) | 1-Hepteno m/z 56 > 42 (≈ 100:55) |
| B11 | Haluros de ácido: ion halocarbonilo X–C≡O⁺ menos estable que R–C≡O⁺ | Carbonilos §2 | Cloruro de butanoilo m/z 71 > 63 |

En las reglas B1 y B6 las etiquetas de átomo (C3, O5…) siempre remiten a los índices del dibujo de la molécula original, aunque el mecanismo ocurra sobre un isómero o un ion intermedio.

### Representación de los mecanismos (v1.3)

Cada paso de la ruta se muestra como un esquema **precursor → ion + neutro**. Sobre el precursor se dibujan:

* **Sitio de ionización**: carga ⊕ y electrón desapareado • sobre el par libre del heteroátomo (ionización n), sobre el enlace π (ionización π) o sobre el enlace σ (ionización σ). En la tarjeta del ion molecular se marca el sitio de menor energía de ionización (n: N > S > I > O > Br > Cl; luego π; luego σ C–C).
* **Flechas curvas** (IUPAC, *curved-arrow notation*): media cabeza (naranja) = desplazamiento de un electrón; completa (azul) = desplazamiento de un par. Origen en el par libre, en el centro del enlace que se rompe o en el H que se transfiere; destino en el átomo o en el punto donde se forma el nuevo enlace.
* **H transferidos** dibujados explícitamente; en los reordenamientos con estado de transición cíclico (McLafferty, McLafferty + 1, eliminaciones de H₂O/HX/ROH/RSH/R–H, ceteno, ion halonio de 5 miembros) la cadena se pliega en el anillo, como en los libros de texto.
* **Enlace que se rompe** sombreado en naranja; índices de átomo de la molécula original en gris.
* Los iones OE⁺• dibujados con su esqueleto neutro se encierran entre corchetes con ⊕•.

Mecanismos con flechas: α-escisión (N, O, S, X; carbonilos; pérdida de H•), escisión inductiva, ion halonio/sulfonio cíclico, escisión alílica, bencílica y propargílica, escisión σ, McLafferty (C=O, C=N, C≡N, C=C, C≡C, arilo) con ambas retenciones de carga y McLafferty + 1, retro-Diels–Alder, eliminaciones de H₂O/H₂S/HX/ROH/RSH/NH₃/R–H, ceteno, •CH₃ de anisoles, α-escisión anular (rutas A, B y C), −CO del ion acilio, reacción del onio y migración del doble enlace (sobre el isómero). Los iones conocidos sólo por su fórmula (−CO de fenoles, −C₂H₂ de aromáticos, apertura de cicloalcanos, −H₂ de cationes) se describen sólo con texto.

Cada ion se clasifica como **OE⁺•** o **EE⁺** (RDB entero/semientero), se verifica la **regla del nitrógeno** y se informan m/z nominal y exacto (restando la masa del electrón).

## 5. Validación

`npm test` compara la predicción con los picos principales de **41 espectros EI de referencia** (valores aproximados, NIST WebBook SRD 69). Además del conjunto original, la v1.1 incorpora ácidos, ésteres, cetonas y aminas cíclicas, un alcohol y un éter cíclicos, pares de isómeros orto/para y cafeína. Se informan tres métricas: picos de referencia predichos (≥ 2 %), acierto del pico base y similitud coseno calculada sobre los m/z de referencia.

| Versión | Conjunto | Picos predichos | Pico base | Coseno medio |
|---|---|---|---|---|
| v1.0 | 41 compuestos | 146/194 (75 %) | 29/41 | — |
| **v1.1** | 41 compuestos | **173/194 (89 %)** | **40/41** | **0,953** |
| v1.1 | 20 originales | 81/91 | 20/20 | 0,960 |
| v1.1 | 21 nuevos | 92/103 | 20/21 | 0,945 |

El único pico base no acertado corresponde al ácido benzoico (predicho 122 > 105; en la referencia 105 > 122, con ambos picos intensos).

**Validación v1.2.** El archivo `test/reference.js` publicado en el repositorio contenía sólo los 20 compuestos originales (los 21 de la v1.1 citados arriba no están en el repositorio). A esos 20 se añadieron 10 compuestos ilustrativos del capítulo de Beauchamp (con datos del propio capítulo y de NIST):

| Conjunto | v1.1 | v1.2 |
|---|---|---|
| 20 originales | 81/91 · base 20/20 | 81/91 · base 20/20 (sin cambios) |
| 10 Beauchamp | 66/85 · base 4/10 | **78/85** · base 4/10 |
| Total 30 | 147/176 | **159/176 (90 %)** · base 24/30 |

| Compuesto | Predicho v1.2 (m/z:%) | Referencia (aprox.) |
|---|---|---|
| 1-Hepteno | 41:100 56:65 55:58 42:43 39:26 69:26 98:20 | 56:100 41:97 55:68 29:56 42:55 70:44 69:31 |
| (E)-2-Hepteno | 55:100 56:70 41:58 69:55 42:54 70:33 98:21 | 55:100 56:90 41:75 69:48 98:44 70:17 |
| 3,4-Dimetilhexano | 57:100 56:60 55:60 41:34 85:34 | 56:100 57:81 43:58 41:43 85:41 |
| 1-Hexanotiol | 47:100 89:31 61:22 118:16 56:15 | 56:100 43:48 41:35 55:35 118:30 |
| 1-Heptino | 39:100 53:40 56:32 43:32 67:29 | 81:100 41:71 55:51 29:46 67:44 |

Picos base aún no acertados: 1-hepteno (41 vs 56, ambos ≈ 100 % en la referencia), 3,4-dimetilhexano (57 vs 56, ambos intensos), 1-butanol (31 vs 56), 1-hexanotiol (47 vs 56), 4-octanona (71 vs 57) y 1-heptino (39 vs 81).

**Compuestos incorporados en la v1.1**

| Compuesto | Predicho v1.1 (m/z:%) | Referencia (aprox.) |
|---|---|---|
| Acetato de etilo | 43:100 61:19 88:19 73:9 | 43:100 45:15 61:14 70:12 88:8 |
| Benzoato de etilo | 105:100 77:64 122:44 150:36 123:23 | 105:100 77:40 122:30 150:25 123:10 |
| Ciclohexanona | 55:100 42:61 41:30 98:22 70:19 | 55:100 42:85 98:35 41:30 69:25 70:15 |
| Ciclopentanona | 55:100 28:58 56:46 84:30 | 55:100 84:50 28:40 56:35 41:30 |
| Ciclohexanol | 57:100 82:30 67:20 100:2 | 57:100 82:45 67:30 44:25 100:3 |
| Piperidina | 84:100 56:69 85:66 57:46 | 84:100 85:60 56:50 57:45 44:35 30:30 |
| Pirrolidina | 43:100 70:78 71:51 42:15 | 43:100 71:45 70:40 42:35 28:30 |
| Tetrahidrofurano | 42:100 41:49 71:26 72:19 | 42:100 41:50 71:45 72:30 27:25 |
| Salicilato de metilo (orto) | **120**:100 121:76 92:73 152:26 | **120**:100 92:60 152:45 121:40 |
| 4-Hidroxibenzoato de metilo (para) | **121**:100 93:54 152:38 | **121**:100 152:50 93:20 |
| 2-Nitrotolueno (orto) | **120**:100 91:71 137:58 92:49 65:44 | **120**:100 65:70 92:55 91:40 137:10 |
| 4-Nitrotolueno (para) | **91**:100 137:83 107:41 | **91**:100 137:70 65:60 107:20 |
| Aspirina | 120:100 92:73 138:66 43:63 180:19 | 120:100 43:45 138:45 92:30 180:10 |
| Cafeína | 194:100 109:86 82:51 55:44 137:33 | 194:100 109:60 55:35 67:35 82:25 |

Los pares orto/para muestran que el motor distingue isómeros de posición mediante el efecto orto.

**Conjunto original (v1.0, sin regresión)**

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
2. No se modelan aún: reordenamientos de esqueleto en policiclos fusionados (la α-escisión en anillos de la v1.1 sólo actúa cuando la doble ruptura separa la molécula, por lo que esteroides, terpenos y decalina siguen sin tratarse bien); retro-Diels–Alder en flavonoides y tetrahidroisoquinolinas; efecto orto en cetonas, amidas y éteres; derivados TMS (m/z 73, 75, 147); pérdidas de CO₂/SO₂/N₂; pérdidas sucesivas de halógenos; iones de carga doble y metaestables.
3. Alquinos (v1.2): el [M−CH₃]⁺ dominante de los 1-alquinos (m/z 81 del 1-heptino) y los iones alquenilo C₃H₅⁺/C₄H₇⁺ requieren reordenamientos de H múltiples no modelados. Tampoco se modelan los iones «alcano radical-catión» de masa par (30, 44, 58…) que Beauchamp atribuye a reordenamientos de alta energía.
4. Algunos picos menores siguen ausentes o subestimados: m/z 31 del terc-butanol, 73 de ácidos alifáticos (escisión γ), 44 y 30 de la piperidina y 67 de la cafeína.
5. Las estructuras de iones de reordenamiento se dibujan con su esqueleto neutro y se rotulan ⁺•; los iones descritos sólo por fórmula (p. ej., C₅H₅⁺) no se dibujan.
6. Los átomos se citan con el índice del dibujo RDKit (base 0).

### Cómo añadir una regla

Cree en `src/ms/rules_ext.js` (o `rules.js`) una función `(g) => eventos[]` que devuelva objetos `{ rule, ruleName, ion: Species, neutral, score, steps, refs, oe }` y agréguela a `EXT_RULES` (o `PRIMARY_RULES`). Las fragmentaciones consecutivas se añaden en `src/ms/secondary.js`. Ejecute `npm test` para verificar que la calibración no se degrada.

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
18. Yates, B. F.; Bouma, W. J.; Radom, L. Distonic radical cations: guidelines for the assessment of their stability. *Tetrahedron* 1986, 42, 6225–6234.
19. Schwarz, H. Some newer aspects of mass spectrometric ortho effects. *Top. Curr. Chem.* 1978, 73, 231–263.
20. Beauchamp, P. *Basics of Mass Spectroscopy* (capítulo del *Spectroscopy Workbook*). California State Polytechnic University, Pomona. Material docente (MS_chapter.pdf).

Licencias de terceros: RDKit (BSD-3), JSME (BSD-3), Express (MIT), Helmet (MIT).

## 7. Historial de versiones

- **v1.3.1** — Corrección: `index.html` se sirve sin caché y `app.js`/`style.css` llevan número de versión (`?v=`), para que el navegador no siga usando la interfaz anterior (que no mostraba los esquemas con flechas) durante la hora de caché.
- **v1.3.0** — Esquemas de mecanismo con flechas curvas (1 e⁻ / 2 e⁻), carga ⊕ y electrón desapareado • en el precursor; sitio de ionización del M⁺•; estados de transición cíclicos plegados; H transferidos explícitos; corchetes ⊕• en iones OE⁺•; leyenda de convenciones en la interfaz. El espectro calculado no cambia.
- **v1.2.0** — Patrones de Beauchamp (MS_chapter.pdf) ausentes en v1.1: migración del doble/triple enlace antes de fragmentarse; escisión propargílica, [M−1]⁺ y McLafferty en alquinos; eliminación de R–H en alcanos ramificados; eliminación de ROH, RSH, NH₃/RNH₂; fragmentación alquénica consecutiva de iones de eliminación; −(H₂O + C₂H₄) en 1-alcanoles; onio cíclicos de 5 y 3 miembros con Cl, Br, I y S; doble McLafferty en cetonas; −H₂ en cationes alilo; ion halocarbonilo. reparto de carga por la regla de Stevenson en el McLafferty de alquenos (1-hepteno 56 > 42). Conjunto de validación ampliado con 10 compuestos del capítulo.
- **v1.1.0** — α-Escisión en anillos con ion distónico (cetonas, alcoholes, aminas, éteres y sulfuros cíclicos); McLafferty + 1 y pérdida de ROH en ésteres; McLafferty desde la cadena alcoxílica con retención de carga en el alqueno; efecto orto (−ROH/−H₂O, −•OH) también en iones secundarios; reacción del onio desde el carbono en oxonios de alcoholes; escisión bencílica con transferencia de H desde N–H/O–H; regla de Stevenson en α-escisiones que expulsan radicales bencílicos; retro-reacción de imidas/ureas cíclicas (cafeína); validación ampliada a 41 compuestos con similitud coseno.
- **v1.0.0** — Versión inicial: 11 reglas primarias, 4 secundarias, 20 compuestos de validación; despliegue PM2/Nginx y Dokploy.
