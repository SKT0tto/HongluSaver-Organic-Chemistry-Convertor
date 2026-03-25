# HongluSaver — Chemistry Formula Converter

**HongluSaver** is a web application for **organic chemistry name–structure interconversion**. It resolves **common names, trade names, and many IUPAC-style strings** to validated structures, then exports multiple representations (2D/3D, SMILES, InChI, molecular formula, properties). A **freehand drawing mode** (Ketcher) turns skeletal formulas into the same pipeline. Display naming is tuned toward **IB Chemistry (International Baccalaureate)** teaching conventions (systematic names preferred over selected retained trivial names where applicable).

**Author:** Henry Lin  
**Repository:** [HongluSaver-Organic-Chemistry-Convertor](https://github.com/SKT0tto/HongluSaver-Organic-Chemistry-Convertor) (example; update if your URL differs)

---

## Table of Contents

1. [Features](#features)  
2. [How It Works (Architecture)](#how-it-works-architecture)  
3. [Prerequisites](#prerequisites)  
4. [Local Installation](#local-installation)  
5. [Running the App](#running-the-app)  
6. [Deployment (Streamlit Community Cloud)](#deployment-streamlit-community-cloud)  
7. [Project Structure](#project-structure)  
8. [Configuration & Caching](#configuration--caching)  
9. [External Services & Data Attribution](#external-services--data-attribution)  
10. [Limitations & Disclaimers](#limitations--disclaimers)  
11. [Development](#development)  
12. [Changelog](#changelog)  
13. [License](#license)  

---

## Features

### Search mode

- **Input:** English compound name (common, trivial, or IUPAC-like string as understood by PubChem).  
- **Output:**  
  - **Preferred display name** — IB-oriented normalization (see [Naming pipeline](#naming-pipeline)).  
  - **PubChem Lexichem / synonym** captions when they differ from the chosen display name.  
  - **2D structure** (RDKit rendering).  
  - **3D structure** — interactive ball-and-stick viewer ([3Dmol.js](https://3dmol.org/)): rotate, zoom, auto-spin.  
  - **Molecular formula** (HTML subscripts), **condensed structural formula** (heuristic for simple acyclic chains; falls back to SMILES when rings/complexity apply).  
  - **SMILES**, **isomeric SMILES** (when distinct), **InChI**, **InChI Key**.  
  - **Molecular properties:** molecular weight, exact mass, heavy atoms, formal charge, H-bond donors/acceptors, rotatable bonds, XLogP (when available from PubChem).  
  - **Synonyms** (short list from PubChem).  
  - **Link** to PubChem compound record (CID) when present.

### Free drawing mode

- Embedded **[Ketcher](https://lifescience.opensource.epam.com/ketcher/)** editor (`streamlit-ketcher`) for skeletal (line-angle) structures.  
- **Validate** drawings with RDKit (connectivity, basic chemistry checks, size limits).  
- **Canonicalize SMILES** and run the **same lookup + display pipeline** as search mode when the structure exists in PubChem; **fallback path** uses RDKit + optional NCI Cactus names when PubChem has no hit.

### Structural isomer explorer

- **Input:** Hill-system style **molecular formula** (e.g. `C4H10O`).  
- **Process:** PubChem PUG REST **formula search** → CID list → property fetch (supports PubChem’s `ConnectivitySMILES` / legacy `CanonicalSMILES`).  
- **Dedup:** canonical SMILES so each **distinct connectivity** appears once (up to a capped batch).  
- **Display:** grid of **2D images**, **IB-oriented names**, molecular weight, SMILES.  
- **Quick-example** formula buttons; results and network steps are **cached** (see [Caching](#configuration--caching)).

### User experience

- **Bilingual UI:** Chinese / English (sidebar toggle).  
- **Sidebar quick examples** for search mode (populate input + optional auto-search pattern).  
- **Branding:** custom CSS, gradient title, benzene-themed logo.  
- **Responsive layout:** wide Streamlit layout; structure + data in columns.

---

## How It Works (Architecture)

### High-level flow

```text
User input (name OR drawn SMILES OR molecular formula)
        │
        ▼
┌───────────────────┐     ┌─────────────────────┐
│  PubChem (names/   │     │  RDKit (parse,       │
│  SMILES, formula)  │     │  2D/3D, descriptors) │
└───────────────────┘     └─────────────────────┘
        │                           │
        └───────────┬───────────────┘
                    ▼
           ┌────────────────┐
           │ Naming layer   │  ← IB table, trivial→systematic map,
           │ (chemistry_     │    regex locant fixes, scoring
           │  engine.py)     │
           └────────────────┘
                    ▼
           Streamlit UI (app.py)
```

### Naming pipeline

HongluSaver does **not** ship a full arbitrary-IUPAC-to-structure parser (that would be OPSIN-class complexity). Instead:

1. **PubChem** resolves the **compound identity** (structure is fixed to a CID / SMILES).  
2. **Display names** are derived from:  
   - A **curated SMILES → IB-friendly name** table for common syllabus molecules.  
   - Otherwise **PubChem Lexichem (IUPAC) name** plus **text normalization**: trivial→systematic dictionary (e.g. *acetic acid* → *ethanoic acid*), **regex** modernization of old locant style (*1-butanol* → *butan-1-ol*).  
3. For **offline-style** fallback (drawn structure not in PubChem), **NCI Cactus** name lines may be fetched and scored with the same **IB-style heuristics**.

This is **pedagogically aligned** with IB expectations; it is **not** a guarantee of the sole official IUPAC name for every edge case worldwide.

### 3D coordinates

3D models are generated in **RDKit** (e.g. ETKDGv3 embedding + force-field optimization when possible), exported as **Mol blocks**, and passed to **3Dmol.js** inside an HTML component. No extra Python package is required for the viewer (CDN script).

### Isomer search

Formula queries use PubChem’s **asynchronous** formula endpoint when needed: initial `202` / `Waiting` + `ListKey`, then **polling** until `IdentifierList` is ready. The implementation retries on transient **SSL/network** errors.

---

## Prerequisites

- **Python** 3.10+ recommended (3.11 tested).  
- **Internet access** for PubChem, Cactus (fallback), and 3Dmol.js CDN.  
- **Git** (optional, for clone/deploy).

On **Linux** servers or some containers, RDKit drawing may need X11-related libraries; for **Streamlit Cloud**, see `packages.txt` in this repo.

---

## Local Installation

### 1. Clone or copy the project

```bash
git clone https://github.com/SKT0tto/HongluSaver-Organic-Chemistry-Convertor.git
cd HongluSaver-Organic-Chemistry-Convertor
```

### 2. Create a virtual environment (recommended)

```bash
python -m venv .venv

# Windows (cmd)
.venv\Scripts\activate

# macOS / Linux
source .venv/bin/activate
```

### 3. Install dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

**Note:** The `requirements.txt` pins **`rdkit`** (official PyPI wheels). Older tutorials may mention `rdkit-pypi`; this project uses the current **`rdkit`** package name.

---

## Running the App

From the project root (with venv activated):

```bash
streamlit run app.py
```

Default URL: **http://localhost:8501**

### Windows shortcut

`HongluSaver.bat` can install dependencies on first run and launch Streamlit (keeps a console window open). Adjust paths if your Python installation is non-standard.

---

## Deployment (Streamlit Community Cloud)

1. Push this repository to **GitHub** (main branch).  
2. On [Streamlit Community Cloud](https://share.streamlit.io/), create a **New app**.  
3. Select the repo, branch, and **main file:** `app.py`.  
4. **Python dependencies** are read from `requirements.txt`.  
5. If RDKit rendering fails on the build image, ensure **`packages.txt`** is present (this repo lists minimal Debian packages for some RDKit backends).  

**Secrets:** The app does not require API keys for PubChem public endpoints. Do not commit secrets; use Streamlit **Secrets** only if you extend the app with private APIs.

**Custom domains / branding:** Follow Streamlit Cloud documentation; author credit appears in-app (footer and sidebar).

---

## Project Structure

```text
.
├── app.py                 # Streamlit UI, i18n, layout, 3D HTML wrapper
├── chemistry_engine.py    # PubChem/RDKit integration, IB naming layers, isomer search
├── requirements.txt       # Python packages
├── packages.txt           # Optional apt packages (e.g. Streamlit Cloud)
├── HongluSaver.bat        # Windows launcher
├── CHANGELOG.md           # Release notes
└── README.md              # This file
```

---

## Configuration & Caching

### Streamlit cache

The app uses `st.cache_data` with a **24-hour TTL** (and entry limits) for:

- Normalized **name lookups**  
- **SMILES-based** compound resolution  
- **2D PNG** generation (by SMILES + dimensions)  
- **3D Mol blocks** (by canonical SMILES)  
- **Isomer search** (by formula string)  

This reduces load on PubChem and speeds up repeat queries. Clear cache from Streamlit’s menu if you need to force fresh data.

### Theme

Light/dark appearance follows the **host environment / Streamlit theme**. Logo text stroke in the sidebar SVG adjusts for dark mode via CSS media queries and Streamlit theme selectors where applicable.

---

## External Services & Data Attribution

| Service | Use |
|--------|-----|
| [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | Compound resolution by name/SMILES/formula; properties, synonyms, systematic names |
| [NCI Cactus](https://cactus.nci.nih.gov/chemical/structure) | Supplementary names when PubChem has no hit (drawing fallback) |
| [3Dmol.js](https://3dmol.org/) | In-browser 3D molecular visualization (loaded from CDN) |
| [Ketcher](https://lifescience.opensource.epam.com/ketcher/) | Molecular editor in drawing mode |

Please respect **PubChem’s usage policies** and cite them in academic work when referring to downloaded data.

---

## Limitations & Disclaimers

- **Not every chemical name** will resolve in PubChem; typos, very new compounds, or proprietary labels may fail.  
- **IUPAC display names** are **heuristic + curated** for teaching; for publication-grade nomenclature, always cross-check primary sources.  
- **Isomer lists** from PubChem are **all entries sharing the formula**, including salts, hydrates, or multi-component records in the database — **chemically filter** results for examination settings if needed.  
- **3D conformers** are computer-generated approximations, not experimental crystal structures.  
- **Condensed formulas** from SMILES are **heuristic** and intentionally simple.  
- **Network**: corporate or school firewalls may block PubChem or CDNs; VPN/proxy configuration may be required for `git push` or API access in some regions.

---

## Development

### Suggested workflow

```bash
pip install -r requirements.txt
streamlit run app.py
```

### Optional quality improvements

- Add **`pytest`** tests for `_ib_normalize`, `_preferred_display_name`, and formula / SMILES edge cases.  
- Pin dependency versions in CI to match Streamlit Cloud.  
- Add **screenshots** to this README (`docs/images/`).

### Code style

The codebase favors **clear separation**: `chemistry_engine.py` for science/network logic, `app.py` for presentation.

---

## Changelog

See **[CHANGELOG.md](./CHANGELOG.md)** for version history.

---

## License

Unless otherwise specified in this repository, you may treat the project as **open for educational use**. If you need a standard license (e.g. MIT), add a `LICENSE` file and update this section accordingly.

---

## Acknowledgements

- **RDKit** — cheminformatics toolkit  
- **PubChem / NIH** — open chemical database  
- **Streamlit** — application framework  
- **EPAM Ketcher & 3Dmol.js teams** — drawing and 3D visualization  

**HongluSaver — Chemistry formula converter** · **Made by Henry Lin**
