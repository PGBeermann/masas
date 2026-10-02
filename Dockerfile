# EI-MS Predictor — imagen para Dokploy (build type: Dockerfile)
FROM node:22-bookworm-slim

ENV NODE_ENV=production \
    PORT=3100 \
    HOST=0.0.0.0 \
    TRUST_PROXY=1 \
    BASE_PATH="" \
    RATE_LIMIT=60

WORKDIR /app

# Dependencias primero (aprovecha la caché de capas de Docker)
COPY package.json package-lock.json ./
RUN npm ci --omit=dev && npm cache clean --force

# Código de la aplicación
COPY --chown=node:node . .

USER node
EXPOSE 3100

# Dokploy/Docker usan este chequeo para marcar el contenedor como sano
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:'+process.env.PORT+(process.env.BASE_PATH||'')+'/api/health').then(r=>r.json()).then(j=>process.exit(j.ok?0:1)).catch(()=>process.exit(1))"

CMD ["node", "server.js"]
