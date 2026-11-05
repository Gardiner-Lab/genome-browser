# Contract: `seed_stock_tags` — Genome Browser ↔ Seed Bank family-tree link

**Contract version: 1**

This document defines the stable, versioned read-side contract between the
**Genome Browser** (which owns and writes seed-stock tags) and the **Seed Bank**
(which reads them to surface a genome-browser link in a line's or plant's
`Family_Tree`). It fulfils the in-repo side of the family-tree integration
(Requirements 11.5 and 11.7) and is the specification the external Seed Bank
source change is implemented against.

Both sides depend only on what is written here. Any change to the row schema,
the key, the ownership rule, or the link shape is a **breaking change** and MUST
bump the contract version above.

---

## 1. Ownership

- The **Genome Browser OWNS and WRITES** the `seed_stock_tags` rows. They live
  in the Genome Browser's own database, `genome.db`, on the `genome-data`
  volume. The Genome Browser inserts a row when a Member adds a `Seed_Stock_Tag`
  to a track or genome (Req 11.2, 11.5) and deletes the row when the Member
  removes the tag (Req 11.7, 11.8). The authoritative writer is
  `genome/tags.py`.
- The **Seed Bank READS** these rows **read-only**, keyed by `stock_id`. It
  never writes, alters, or deletes them.
- **Removing a row removes the family-tree link.** When the Genome Browser
  deletes a `seed_stock_tags` row (Member removes the tag), the corresponding
  link disappears from the Seed Bank's `Family_Tree` on its next read, because
  the Seed Bank's view is derived from these rows and holds no copy of its own
  (Req 11.7).

---

## 2. Database location

| Property | Value |
|----------|-------|
| Volume | `genome-data` |
| Mounted into the `seeds` container at | `/genome` **read-only** (`genome-data:/genome:ro`) |
| Database file | `genome.db` (i.e. `/genome/genome.db` from inside the Seed Bank container) |
| Table | `seed_stock_tags` |

The read-only mount is declared in `docker-compose.yml` on the `seeds` service
(done as part of task 16.4). The Seed Bank opens `/genome/genome.db` read-only;
it must tolerate the file being absent (the Genome Browser creates it on first
startup) and treat "absent" as "no family-tree links from the Genome Browser".

---

## 3. Row schema

Source of truth: `genome/schema.sql`. The table is created idempotently by
`genome/db.py` at startup.

```sql
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
```

| Column | Type | Meaning |
|--------|------|---------|
| `id` | TEXT (UUID hex) | Primary key of the tag row. Opaque to the Seed Bank. |
| `target_type` | TEXT | What the tag is attached to: `'track'` or `'genome'`. Tells the Seed Bank which kind of genome-browser record to link to (see §5). |
| `target_id` | TEXT | The id of that track or reference genome in the Genome Browser. This is the value the link path is built from. |
| `stock_kind` | TEXT | The Seed Bank record kind: `'line'`, `'lot'` or `'plant'`. Decides whether the row is a `Family_Tree` node (see §4). |
| `stock_id` | TEXT | **The key.** The external Seed Bank id of the line / lot / plant. For a line or plant this id *is* the `Family_Tree` key the Seed Bank joins on. |
| `verified_at` | TEXT | Timestamp the tag was recorded (`datetime('now')`). Informational; the Seed Bank may display it but must not depend on its format beyond "SQLite datetime text". |
| `created_by` | TEXT | The Member who added the tag. Informational. |

**Uniqueness:** `UNIQUE (target_type, target_id, stock_id)` — a given stock may
be tagged on a given target at most once (Req 11.4). The Seed Bank may therefore
see more than one row per `stock_id` (the same line tagged on several
tracks/genomes); each such row is a distinct family-tree link.

### Reading by key

```sql
SELECT target_type, target_id, stock_kind, stock_id, verified_at, created_by
FROM seed_stock_tags
WHERE stock_id = ? AND stock_kind IN ('line', 'plant');
```

---

## 4. Where it surfaces in the Seed Bank

- A row with `stock_kind = 'line'` or `stock_kind = 'plant'` surfaces as a link
  in the **`Family_Tree` of the line or plant identified by `stock_id`**
  (Req 11.5). These are the only two kinds that are `Family_Tree` nodes.
- A row with `stock_kind = 'lot'` is **recorded** by the Genome Browser (a seed
  lot can be tagged), but a **seed lot is not a `Family_Tree` node**. The Seed
  Bank should not expect to place a `'lot'` row in the family tree; it may
  surface it elsewhere (e.g. on the seed lot's own page) but that is outside
  this contract. The Seed Bank's family-tree read SHOULD filter to
  `stock_kind IN ('line','plant')`.

---

## 5. Resolving a tag to a navigable Genome Browser link

The Genome Browser is mounted at **`/genome`** in the suite (its `URL_PREFIX`).
The Seed Bank reaches the Genome Browser through the shared suite URL space
(the same convention the Genome Browser uses for `SEEDS_URL=/seeds` when it
links *out* to the Seed Bank, see `genome/tags.py`). From a `seed_stock_tags`
row the Seed Bank builds the browser link from `target_type` and `target_id`:

| `target_type` | Link shape | Links to |
|---------------|-----------|----------|
| `'track'` | `/genome/tracks/{target_id}` | the tagged track in the Genome Browser |
| `'genome'` | `/genome/genomes/{target_id}` | the tagged reference genome |

Notes to match the Genome Browser exactly:

- `/genome` is the Genome Browser mount. If the Seed Bank is configured with a
  genome-browser base URL, it is that base; otherwise default to `/genome`.
- `target_id` MUST be URL-encoded as a single path segment, mirroring the
  Genome Browser's own `_join()` helper in `genome/tags.py`, which quotes each
  path part with `urllib.parse.quote(part, safe='')`. This keeps ids containing
  spaces or slashes forming a valid link.
- The link is read-only navigation; the Seed Bank never calls a Genome Browser
  write endpoint from the family tree.

For symmetry: the *outbound* direction (Genome Browser → Seed Bank) already
exists in `genome/tags.py` `_seed_stock_link()`, which builds
`{SEEDS_URL}/{lines|lots|plants}/{stock_id}`. This contract defines the inbound
direction (Seed Bank → Genome Browser) that the family tree needs.

---

## 6. EXTERNAL HANDOFF — Seed Bank source change (follow-up deliverable)

The Seed Bank ships as a **prebuilt image** (`ghcr.io/gardiner-lab/seed-bank`)
whose source is **not in this repository**. The actual code that:

1. opens `/genome/genome.db` read-only,
2. reads `seed_stock_tags` by `stock_id` (filtered to `stock_kind IN
   ('line','plant')`), and
3. renders the navigable Genome Browser link (per §5) in the `Family_Tree` of
   the matching line or plant,

is therefore an **EXTERNAL follow-up deliverable**, implemented in the Seed Bank
codebase against **this contract (version 1)**.

### In-repo deliverables (complete)

- ✅ **Compose mount** — `genome-data:/genome:ro` added to the `seeds` service in
  `docker-compose.yml`, giving the Seed Bank read-only access to the Genome
  Browser's `genome.db`.
- ✅ **This contract** — `genome/CONTRACT_seed_stock_tags.md`, the stable,
  versioned schema and link convention both sides depend on.
- ✅ **Write side** — `genome/tags.py` and `genome/schema.sql` already own and
  write the `seed_stock_tags` rows this contract exposes.

### External side (not in this repo)

- ⛔ **Seed Bank code** — reading `seed_stock_tags` and rendering the family-tree
  link. Tracked as a follow-up against the Seed Bank image; it consumes the
  mount and the schema defined here without any further change on the Genome
  Browser side, as long as the contract version is unchanged.
