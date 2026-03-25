# Changelog

All notable changes to **HongluSaver - Chemistry Formula Converter** will be documented in this file.

---

## [1.2.1] - 2026-03-26

### Added
- **Special thanks** to Henry Lu (aka Honglu) in README, app footer, and GitHub Pages landing page, with links to the Streamlit app and GitHub repository.

---

## [1.2.0] - 2026-03-25

### Added
- **3D Molecular Viewer**: Interactive 3D ball-and-stick model powered by 3Dmol.js — drag to rotate, scroll to zoom, auto-spin. Available as a tab alongside the 2D structure.
- **Structural Isomer Generator**: New "Isomers" tab — enter a molecular formula (e.g. C₄H₁₀O) to find all structural isomers from PubChem, displayed in a 3-column grid with 2D structures and IB-compliant names.
- 7 quick-access buttons for common IB molecular formulas.
- Cached isomer search results (24 h TTL).

### Fixed
- PubChem API now returns `ConnectivitySMILES` instead of `CanonicalSMILES`; updated field reading to support both.
- Isomer example buttons no longer crash with `StreamlitAPIException` (switched to `on_click` callback).
- Added network retry mechanism (`_pug_get`) for PubChem polling to handle transient SSL/connection drops.

---

## [1.1.1] - 2026-03-25

### Improved
- Sidebar logo: regular hexagon geometry (equal side lengths), SVG-native `<text>` for crisp rendering.
- Dual-layer SVG text with white stroke (light mode) / dark stroke (dark mode) for maximum readability.

---

## [1.1.0] - 2026-03-25

### Added
- Bilingual UI: Chinese / English toggle in sidebar.
- Free Drawing Mode: Ketcher molecular editor for sketching skeletal formulas → auto IUPAC name generation.
- Sidebar branding: benzene ring SVG logo with "HongluSaver" text.
- Quick example buttons in sidebar for common compounds.
- Result caching (`st.cache_data`, 24 h TTL) for compound lookups, SMILES lookups, and 2D image rendering.

### Improved
- Animated gradient title with Orbitron font.
- IB Chemistry compliant naming:
  - Hard-coded SMILES → IB name dictionary (~120 compounds).
  - Trivial-to-systematic mapping (e.g. acetic acid → ethanoic acid).
  - Regex locant correction (1-butanol → butan-1-ol).
  - IB-aware scoring system for name ranking.
- Removed OPSIN API dependency for faster queries.

---

## [1.0.0] - 2026-03-25

### Added
- Initial release: IUPAC name → structure converter.
- Compound lookup via PubChem API.
- Output: molecular formula, condensed structural formula, SMILES, InChI, InChI Key, 2D structure image.
- Molecular properties: MW, exact mass, H-bond donors/acceptors, rotatable bonds, XLogP.
- Synonyms display from PubChem.
- Streamlit web UI with responsive layout.
