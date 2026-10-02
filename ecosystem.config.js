// PM2 — pm2 start ecosystem.config.js && pm2 save
module.exports = {
  apps: [{
    name: 'eims-predictor',
    script: 'server.js',
    instances: 2,              // RDKit-WASM es síncrono: varias instancias reparten la carga
    exec_mode: 'cluster',
    max_memory_restart: '400M',
    env: { NODE_ENV: 'production', PORT: 3100, HOST: '127.0.0.1', BASE_PATH: '', RATE_LIMIT: 60 }
  }]
};
