# HongluSaver — Chemistry Formula Converter

**HongluSaver** is a web application for **organic chemistry name–structure interconversion**. It resolves **common names, trade names, and many IUPAC-style strings** to validated structures, then exports multiple representations (2D/3D, SMILES, InChI, molecular formula, properties). A **freehand drawing mode** (Ketcher) turns skeletal formulas into the same pipeline. Display naming is tuned toward **IB Chemistry (International Baccalaureate)** teaching conventions (systematic names preferred over selected retained trivial names where applicable).

**Author:** Henry Lin  
**Live app:** [https://honglusaver.streamlit.app/](https://honglusaver.streamlit.app/)  
**Repository:** [github.com/SKT0tto/HongluSaver-Organic-Chemistry-Converter](https://github.com/SKT0tto/HongluSaver-Organic-Chemistry-Converter)

---

## Table of Contents

1. [Live demo](#live-demo)  
2. [Features](#features)  
3. [How It Works (Architecture)](#how-it-works-architecture)  
4. [Prerequisites](#prerequisites)  
5. [Local Installation](#local-installation)  
6. [Running the App](#running-the-app)  
7. [Deployment (Streamlit Community Cloud)](#deployment-streamlit-community-cloud)  
8. [Project Structure](#project-structure)  
9. [Configuration & Caching](#configuration--caching)  
10. [External Services & Data Attribution](#external-services--data-attribution)  
11. [Limitations & Disclaimers](#limitations--disclaimers)  
12. [Development](#development)  
13. [Changelog](#changelog)  
14. [License](#license)  

---

## Live demo

The production build is hosted on **Streamlit Community Cloud**:

**[https://honglusaver.streamlit.app/](https://honglusaver.streamlit.app/)**

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

`HongluSaver.bat` installs dependencies on first run and launches Streamlit (keeps a console window open).

---

## Deployment (Streamlit Community Cloud)

The public instance runs at **[honglusaver.streamlit.app](https://honglusaver.streamlit.app/)**, deployed from this GitHub repository via [Streamlit Community Cloud](https://share.streamlit.io/).

**Stack notes for redeployment or forks:**

- Main entrypoint: **`app.py`**.  
- Python packages: **`requirements.txt`**.  
- **`packages.txt`** lists optional Debian packages used on the hosted image for RDKit-related rendering.  

The app uses only **public** PubChem endpoints; no API keys are required for the current feature set.

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

This reduces load on PubChem and speeds up repeat queries. Cached entries can be cleared from the Streamlit app menu when a full refresh is desired.

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

```bash
pip install -r requirements.txt
streamlit run app.py
```

Layout convention: **`chemistry_engine.py`** holds PubChem/RDKit integration and naming logic; **`app.py`** holds the Streamlit UI and presentation.

---

## Changelog

See **[CHANGELOG.md](./CHANGELOG.md)** for version history.

---

## License

Copyright © Henry Lin. All rights reserved. The application is provided for **educational use**.

---

## Acknowledgements

- **RDKit** — cheminformatics toolkit  
- **PubChem / NIH** — open chemical database  
- **Streamlit** — application framework  
- **EPAM Ketcher & 3Dmol.js teams** — drawing and 3D visualization  

**HongluSaver — Chemistry formula converter** · **Made by Henry Lin**
