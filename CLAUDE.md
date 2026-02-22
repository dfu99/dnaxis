# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DNAxiS is a browser-based web application for DNA origami design of shapes with axial symmetry, built with Python 3.8 and Flask. Hosted at caddna.cs.duke.edu.

## Running Locally

```bash
# Prerequisites: Redis must be running
redis-server

# Setup
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# Run
python DNAxiS-websrv.py
# Access at http://127.0.0.1:5000
```

There is no test suite, linter configuration, or CI/CD pipeline.

## Architecture

### Request Flow

The app follows a multi-step wizard workflow with a progress bar (wizard_step 1-5):

1. **Home** (`/`) - Landing page
2. **STL Upload** (`/stl`) - Optional 3D shape upload
3. **Submission** (`/submission`) - Ring drawing canvas, scaffold selection, advanced options (wizard_step=1)
4. **Connections** (`/connections`) - Connection graph editor (wizard_step=2)
5. **Pathway** (`/pathway`) - Routing pathway editor (wizard_step=3)
6. **Verification** (`/pre-submit`) - Pre-flight checks (wizard_step=4)
7. **Processing** (`/process`) - Runs DNA origami generation via `driver.py` (wizard_step=5)
8. **Results** (`/results`) - Download output files (sequences CSV, topology, configuration) (wizard_step=5)

Data passes between steps via Flask-Session backed by Redis. The `require_session_keys()` decorator in `app/routes/session_helpers.py` redirects users to home if required session data is missing.

### Code Organization

- **`app/__init__.py`** - App factory (`create_app()`) configuring Flask-Mail, Flask-Session, and registering the Blueprint. Exports module-level `app` for backward compatibility.
- **`app/scaffold_library.py`** - Single source of truth for `SCAFFOLD_LENGTHS` dict and `ALLOWED_EXTENSIONS` set.
- **`app/routes/`** - Blueprint-based route handlers (merged from former `views/` and `forms/` packages):
  - `__init__.py` - Blueprint definition, imports all route modules
  - `info.py` - `/`, `/test`, `/tutorial`
  - `stl.py` - `/stl`, `/show-example`, `/upload-stl`, `/stl2nodes`
  - `submission.py` - `/submission`, `/submit`, `/update_sequences`
  - `connections.py` - `/connections`, `/upload_connections`
  - `pathway.py` - `/pathway`, `/upload_pathway`
  - `verify.py` - `/pre-submit`, `/verified`
  - `results.py` - `/processing`, `/results`, `/download`, `/download/<path>`
  - `process.py` - `/process` (heavy computation route)
  - `utils.py` - Graph helpers (connectivity, cycle detection, edge validation)
  - `session_helpers.py` - `require_session_keys()` decorator
- **`app/routing/`** - Core DNA routing algorithms (~9,500 LOC). Key modules:
  - `sequence.py` - Scaffold sequence selection
  - `sym_origami.py` - Main origami data structure
  - `crossover.py` - Crossover placement
  - `noncyclic.py`, `modules.py` - Routing strategies
  - `helper/` - Math, trig, rotation matrices, strand navigation utilities
  - `filehandler/` - I/O for ring, pathway, connection, and transaction files
  - `mode_symmetric/` and `mode_asymmetric/` - Routing modes
  - `score/` - Scoring algorithms for routing decisions
- **`app/sequences/`** - Scaffold sequence library (text files + `seqlib.json` metadata). `update.py` handles sequence refresh logic.
- **`app/templates/`** - Jinja2 HTML templates using Bootstrap 5.3 via CDN
- **`app/static/`** - JS, CSS, example STL files
  - `js/draw.js`, `js/cross.js`, `js/preview.js` - Canvas drawing code (DO NOT MODIFY)
  - `js/common.js` - Shared utilities (setError, collapsible toggle)
  - `js/submission.js` - Scaffold selection logic
  - `js/stl-upload.js` - STL file upload/preview
  - `js/connections-page.js`, `js/pathway-page.js`, `js/verify-page.js` - Page-specific AJAX handlers
  - `styles/style.css` - Custom styles (Bootstrap handles layout/nav/footer)

### Configuration

`config.ini` is the central configuration file, loaded globally via `app/config.py` into module-level variables. Key sections: `GENERAL` (available sequences), `NICKING`, `CROSSOVER`, `HELICAL GEOMETRY`, `SIMULATED ANNEALING`, `DIRECTORY`.

Adding a new scaffold sequence requires:
1. Adding the sequence text file to `app/sequences/`
2. Adding the name to `AVAIL_SEQUENCES` in `config.ini`
3. Adding length to `SCAFFOLD_LENGTHS` dict in `app/scaffold_library.py`
4. Updating `app/sequences/seqlib.json`

### Key Patterns

- **App factory + Blueprint**: `create_app()` in `app/__init__.py` creates the Flask app. All routes use `@bp.route` on a single Blueprint registered in the factory.
- **Session-driven state**: All inter-step data stored in Flask session (ring data, connections, pathway, scaffold selections)
- **Job output**: Written to `jobs/<timestamp>/` directories (CSV, TOP, CONF files, packaged as ZIP)
- **Frontend**: Bootstrap 5.3 for layout/responsive design, canvas drawing via jQuery (draw.js, cross.js), 3D preview via Three.js, AJAX for scaffold selection and step submissions
- **Template variables**: Template-specific data (e.g. `circleCoords`) passed via `<script>` blocks setting globals, then external JS reads them
