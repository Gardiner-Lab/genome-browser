-- Genome Browser data model (SQLite, WAL mode; applied by genome/db.py init_db).
--
-- All ids are UUID text, matching the suite. Labs and chair groups are
-- referenced by their suite.db ids/names (read from labs.py), never copied.
-- Every statement is idempotent (IF NOT EXISTS) so startup schema-init is safe
-- to run repeatedly against an existing database.

-- Reference genomes ------------------------------------------------------
CREATE TABLE IF NOT EXISTS reference_genomes (
  id            TEXT PRIMARY KEY,          -- unique, immutable (Req 2.8)
  name          TEXT NOT NULL,
  assembly_id   TEXT NOT NULL,
  owner_lab     TEXT NOT NULL,             -- suite.db lab id
  genome_visibility TEXT NOT NULL,         -- 'chair_group_default' | 'owner_lab_only'
  created_by    TEXT NOT NULL,
  created_at    TEXT NOT NULL DEFAULT (datetime('now')),
  UNIQUE (owner_lab, assembly_id)          -- duplicate assembly id per lab (Req 2.6)
);

CREATE TABLE IF NOT EXISTS chrom_sizes (
  genome_id     TEXT NOT NULL REFERENCES reference_genomes(id) ON DELETE CASCADE,
  chrom         TEXT NOT NULL,             -- exact, case-sensitive name
  length_bp     INTEGER NOT NULL,
  PRIMARY KEY (genome_id, chrom)
);

-- Genome sharing: chair groups ONLY (Req 12.3, 12.4) -------------------
CREATE TABLE IF NOT EXISTS genome_shares (
  genome_id       TEXT NOT NULL REFERENCES reference_genomes(id) ON DELETE CASCADE,
  chair_group_id  TEXT NOT NULL,           -- suite.db chair_groups id
  granted_by      TEXT NOT NULL,
  granted_at      TEXT NOT NULL DEFAULT (datetime('now')),
  PRIMARY KEY (genome_id, chair_group_id)
);

-- Tracks -----------------------------------------------------------------
CREATE TABLE IF NOT EXISTS tracks (
  id            TEXT PRIMARY KEY,
  genome_id     TEXT NOT NULL REFERENCES reference_genomes(id) ON DELETE CASCADE,
  name          TEXT NOT NULL,
  type          TEXT NOT NULL,             -- 'bigwig' | 'bed' | 'gff'
  owner_lab     TEXT NOT NULL,
  visibility    TEXT NOT NULL DEFAULT 'lab_only',  -- 'lab_only' | 'labs' | 'chair_group'
  track_status  TEXT NOT NULL DEFAULT 'ready',     -- 'converting' | 'ready' | 'failed'
  error         TEXT,                      -- failure reason when track_status='failed'
  file_path     TEXT,                      -- stored BigWig/feature file; NULL until ready
  created_by    TEXT NOT NULL,
  created_at    TEXT NOT NULL DEFAULT (datetime('now'))
);

-- Per-track sharing: specific labs OR the owner lab's chair group (Req 7)
CREATE TABLE IF NOT EXISTS track_shares (
  track_id        TEXT NOT NULL REFERENCES tracks(id) ON DELETE CASCADE,
  target_type     TEXT NOT NULL,           -- 'lab' | 'chair_group'
  target_id       TEXT NOT NULL,           -- suite.db lab id or chair_group id
  granted_by      TEXT NOT NULL,
  granted_at      TEXT NOT NULL DEFAULT (datetime('now')),
  PRIMARY KEY (track_id, target_type, target_id)
);

-- Durable conversion queue (survives restarts) (Req 13.12, 13.18, 13.20)
CREATE TABLE IF NOT EXISTS conversion_jobs (
  id            TEXT PRIMARY KEY,
  track_id      TEXT NOT NULL REFERENCES tracks(id) ON DELETE CASCADE,
  kind          TEXT NOT NULL,             -- 'bam' | 'wig'
  input_path    TEXT NOT NULL,             -- scratch original
  state         TEXT NOT NULL DEFAULT 'queued',  -- 'queued' | 'running' | 'done' | 'failed'
  attempts      INTEGER NOT NULL DEFAULT 0,
  error         TEXT,
  created_at    TEXT NOT NULL DEFAULT (datetime('now')),
  started_at    TEXT,
  finished_at   TEXT
);
CREATE INDEX IF NOT EXISTS conversion_jobs_queue ON conversion_jobs (state, created_at);

-- Browser view (persisted per member + genome) (Req 6, 8) ---------------
CREATE TABLE IF NOT EXISTS browser_views (
  username      TEXT NOT NULL,
  genome_id     TEXT NOT NULL REFERENCES reference_genomes(id) ON DELETE CASCADE,
  locus         TEXT,                      -- 'Chr1:1000-5000'
  enabled_json  TEXT NOT NULL DEFAULT '[]',-- ordered track ids
  options_json  TEXT NOT NULL DEFAULT '{}',-- per-track colour/height/signal range
  updated_at    TEXT NOT NULL DEFAULT (datetime('now')),
  PRIMARY KEY (username, genome_id)
);

-- Plasmid / seed-stock tags (Req 11) ------------------------------------
CREATE TABLE IF NOT EXISTS plasmid_tags (
  id            TEXT PRIMARY KEY,
  target_type   TEXT NOT NULL,             -- 'track' | 'genome'
  target_id     TEXT NOT NULL,
  plasmid_id    TEXT NOT NULL,             -- external MoClo/Editor id
  verified_at   TEXT NOT NULL DEFAULT (datetime('now')),
  created_by    TEXT NOT NULL,
  UNIQUE (target_type, target_id, plasmid_id)   -- duplicate tag (Req 11.4)
);

CREATE TABLE IF NOT EXISTS seed_stock_tags (
  id            TEXT PRIMARY KEY,
  target_type   TEXT NOT NULL,             -- 'track' | 'genome'
  target_id     TEXT NOT NULL,
  stock_kind    TEXT NOT NULL,             -- 'line' | 'lot' | 'plant'
  stock_id      TEXT NOT NULL,             -- external Seed Bank id (family-tree key)
  verified_at   TEXT NOT NULL DEFAULT (datetime('now')),
  created_by    TEXT NOT NULL,
  UNIQUE (target_type, target_id, stock_id)
);
