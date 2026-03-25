"""
Core chemistry engine: resolves IUPAC names to various chemical representations
using PubChem and RDKit. Preferred display names favor IUPAC substitutive/systematic
forms when identifiable, ranked from PubChem/Cactus synonyms (no per-name web parsing).
"""

import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import pubchempy as pcp
from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, AllChem

_HTTP_UA = "HongluSaver/1.0 (educational; +https://pubchem.ncbi.nlm.nih.gov/)"


@dataclass
class CompoundInfo:
    iupac_name: str
    original_query: str = ""
    pubchem_lexichem_name: str = ""
    common_name: str = ""
    molecular_formula: str = ""
    molecular_weight: float = 0.0
    canonical_smiles: str = ""
    isomeric_smiles: str = ""
    inchi: str = ""
    inchi_key: str = ""
    exact_mass: float = 0.0
    xlogp: Optional[float] = None
    hbond_donor_count: int = 0
    hbond_acceptor_count: int = 0
    rotatable_bond_count: int = 0
    heavy_atom_count: int = 0
    charge: int = 0
    synonyms: list = field(default_factory=list)
    cid: int = 0


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  LAYER 1 — IB Chemistry hard-coded SMILES → IB-standard name table     ║
# ║  Covers all functional groups in IB Chem SL/HL syllabus.               ║
# ║  Keys are RDKit canonical SMILES (built at import time).               ║
# ╚══════════════════════════════════════════════════════════════════════════╝

_IB_TABLE_RAW: dict[str, str] = {
    # ── Alkanes ──
    "C": "methane",
    "CC": "ethane",
    "CCC": "propane",
    "CCCC": "butane",
    "CC(C)C": "2-methylpropane",
    "CCCCC": "pentane",
    "CC(C)CC": "2-methylbutane",
    "CC(C)(C)C": "2,2-dimethylpropane",
    "CCCCCC": "hexane",
    "CC(C)CCC": "2-methylpentane",
    "CCC(C)CC": "3-methylpentane",
    "CC(C)(C)CC": "2,2-dimethylbutane",
    "CC(C)C(C)C": "2,3-dimethylbutane",
    "CCCCCCC": "heptane",
    "CCCCCCCC": "octane",
    # ── Cycloalkanes ──
    "C1CC1": "cyclopropane",
    "C1CCC1": "cyclobutane",
    "C1CCCC1": "cyclopentane",
    "C1CCCCC1": "cyclohexane",
    # ── Alkenes ──
    "C=C": "ethene",
    "CC=C": "propene",
    "C=CCC": "but-1-ene",
    "CC=CC": "but-2-ene",
    "CC(=C)C": "2-methylpropene",
    "C=CCCC": "pent-1-ene",
    "CC=CCC": "pent-2-ene",
    "C=CCCCC": "hex-1-ene",
    # ── Alkynes ──
    "C#C": "ethyne",
    "CC#C": "propyne",
    "C#CCC": "but-1-yne",
    "CC#CC": "but-2-yne",
    "C#CCCC": "pent-1-yne",
    "CC#CCC": "pent-2-yne",
    # ── Alcohols ──
    "CO": "methanol",
    "CCO": "ethanol",
    "CCCO": "propan-1-ol",
    "CC(O)C": "propan-2-ol",
    "CCCCO": "butan-1-ol",
    "CCC(O)C": "butan-2-ol",
    "CC(C)CO": "2-methylpropan-1-ol",
    "CC(C)(O)C": "2-methylpropan-2-ol",
    "CCCCCO": "pentan-1-ol",
    "CCCC(O)C": "pentan-2-ol",
    "CCC(O)CC": "pentan-3-ol",
    "CCCCCCO": "hexan-1-ol",
    # ── Aldehydes ──
    "C=O": "methanal",
    "CC=O": "ethanal",
    "CCC=O": "propanal",
    "CCCC=O": "butanal",
    "CCCCC=O": "pentanal",
    "CCCCCC=O": "hexanal",
    # ── Ketones ──
    "CC(C)=O": "propanone",
    "CCC(C)=O": "butanone",
    "CCCC(C)=O": "pentan-2-one",
    "CCC(CC)=O": "pentan-3-one",
    "CCCCC(C)=O": "hexan-2-one",
    "CCCC(CC)=O": "hexan-3-one",
    # ── Carboxylic acids ──
    "OC=O": "methanoic acid",
    "CC(O)=O": "ethanoic acid",
    "CCC(O)=O": "propanoic acid",
    "CCCC(O)=O": "butanoic acid",
    "CCCCC(O)=O": "pentanoic acid",
    "CCCCCC(O)=O": "hexanoic acid",
    # ── Esters (methyl / ethyl of C1-C4 acids) ──
    "COC=O": "methyl methanoate",
    "COC(C)=O": "methyl ethanoate",
    "COC(=O)CC": "methyl propanoate",
    "COC(=O)CCC": "methyl butanoate",
    "CCOC=O": "ethyl methanoate",
    "CCOC(C)=O": "ethyl ethanoate",
    "CCOC(=O)CC": "ethyl propanoate",
    "CCOC(=O)CCC": "ethyl butanoate",
    # ── Amines ──
    "CN": "methanamine",
    "CCN": "ethanamine",
    "CCCN": "propan-1-amine",
    "CC(N)C": "propan-2-amine",
    "CCCCN": "butan-1-amine",
    "CCC(N)C": "butan-2-amine",
    # ── Amides ──
    "NC=O": "methanamide",
    "CC(N)=O": "ethanamide",
    "CCC(N)=O": "propanamide",
    "CCCC(N)=O": "butanamide",
    # ── Halogenoalkanes ──
    "CCl": "chloromethane",
    "CBr": "bromomethane",
    "CI": "iodomethane",
    "CF": "fluoromethane",
    "ClCC": "chloroethane",
    "BrCC": "bromoethane",
    "ClCCC": "1-chloropropane",
    "CC(Cl)C": "2-chloropropane",
    "BrCCC": "1-bromopropane",
    "CC(Br)C": "2-bromopropane",
    "ClCCCC": "1-chlorobutane",
    "CC(Cl)CC": "2-chlorobutane",
    "BrCCCC": "1-bromobutane",
    "CC(Br)CC": "2-bromobutane",
    # ── Benzene ring compounds (IB-accepted retained + systematic) ──
    "c1ccccc1": "benzene",
    "Cc1ccccc1": "methylbenzene",
    "CCc1ccccc1": "ethylbenzene",
    "Oc1ccccc1": "phenol",
    "Nc1ccccc1": "phenylamine",
    "OC(=O)c1ccccc1": "benzoic acid",
    "Clc1ccccc1": "chlorobenzene",
    "Brc1ccccc1": "bromobenzene",
    "Oc1ccc(O)cc1": "benzene-1,4-diol",
    # ── Amino acids (IB HL) ──
    "NCC(O)=O": "2-aminoethanoic acid",
    "NC(C)C(O)=O": "2-aminopropanoic acid",
    # ── Diols / dicarboxylic (common IB) ──
    "OCC(O)CO": "propane-1,2,3-triol",
    "OC(=O)C(O)=O": "ethanedioic acid",
    "OC(=O)CC(O)=O": "propanedioic acid",
    "OC(=O)CCC(O)=O": "butanedioic acid",
}

_IB_TABLE: dict[str, str] = {}


def _build_ib_table():
    for raw_smi, ib_name in _IB_TABLE_RAW.items():
        mol = Chem.MolFromSmiles(raw_smi)
        if mol is not None:
            canon = Chem.MolToSmiles(mol, canonical=True)
            _IB_TABLE[canon] = ib_name


_build_ib_table()


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  LAYER 2 — Trivial / retained name → IB substitutive name mapping      ║
# ║  Applied as text replacement on any candidate name string.              ║
# ╚══════════════════════════════════════════════════════════════════════════╝

_TRIVIAL_TO_IB: dict[str, str] = {
    "formic acid": "methanoic acid",
    "acetic acid": "ethanoic acid",
    "propionic acid": "propanoic acid",
    "butyric acid": "butanoic acid",
    "valeric acid": "pentanoic acid",
    "caproic acid": "hexanoic acid",
    "formaldehyde": "methanal",
    "acetaldehyde": "ethanal",
    "acetone": "propanone",
    "ethylene": "ethene",
    "propylene": "propene",
    "acetylene": "ethyne",
    "toluene": "methylbenzene",
    "aniline": "phenylamine",
    "glycerol": "propane-1,2,3-triol",
    "glycerine": "propane-1,2,3-triol",
    "glycine": "2-aminoethanoic acid",
    "alanine": "2-aminopropanoic acid",
    "oxalic acid": "ethanedioic acid",
    "malonic acid": "propanedioic acid",
    "succinic acid": "butanedioic acid",
    "methyl alcohol": "methanol",
    "ethyl alcohol": "ethanol",
    "propyl alcohol": "propan-1-ol",
    "isopropyl alcohol": "propan-2-ol",
    "butyl alcohol": "butan-1-ol",
    "sec-butyl alcohol": "butan-2-ol",
    "isobutyl alcohol": "2-methylpropan-1-ol",
    "tert-butyl alcohol": "2-methylpropan-2-ol",
    "methyl acetate": "methyl ethanoate",
    "ethyl acetate": "ethyl ethanoate",
    "methyl formate": "methyl methanoate",
    "ethyl formate": "ethyl methanoate",
    "formamide": "methanamide",
    "acetamide": "ethanamide",
    "methylamine": "methanamine",
    "ethylamine": "ethanamine",
    "propylamine": "propan-1-amine",
    "isopropylamine": "propan-2-amine",
    "butylamine": "butan-1-amine",
    "methyl chloride": "chloromethane",
    "methyl bromide": "bromomethane",
    "ethyl chloride": "chloroethane",
    "ethyl bromide": "bromoethane",
    "vinyl chloride": "chloroethene",
    "chloroform": "trichloromethane",
    "iodoform": "triiodomethane",
    "bromoform": "tribromomethane",
}


def _apply_trivial_mapping(name: str) -> str:
    """If the whole name (case-insensitive) is a known trivial name, return IB form."""
    return _TRIVIAL_TO_IB.get(name.lower().strip(), name)


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  LAYER 3 — Regex locant & suffix corrections (IB / 2013 IUPAC style)   ║
# ║  Converts old-style "1-butanol" → "butan-1-ol" etc.                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝

_PARENT_STEMS = {
    "methan": "methan", "ethan": "ethan", "propan": "propan", "butan": "butan",
    "pentan": "pentan", "hexan": "hexan", "heptan": "heptan", "octan": "octan",
    "nonan": "nonan", "decan": "decan",
}

_SUFFIX_MAP = {
    "ol": "ol", "one": "one", "amine": "amine", "amide": "amide", "thiol": "thiol",
}


def _fix_old_locant_suffix(name: str) -> str:
    """
    Convert old-style locant placement to IB/2013 style:
      "1-butanol"   → "butan-1-ol"
      "2-pentanone"  → "pentan-2-one"
      "1-propanol"   → "propan-1-ol"
      "2-butanamine" → "butan-2-amine"
    Only applies to simple cases where it's unambiguous.
    """
    n = name.strip()

    # Pattern: (locants)-(parent+suffix)  e.g. "2-pentanone", "1,3-butanediol"
    for suffix_long, suffix_short in [
        ("anol", "ol"), ("anone", "one"), ("anamine", "amine"),
        ("anamide", "amide"), ("anthiol", "thiol"),
    ]:
        m = re.match(
            rf"^(\d+(?:,\d+)*)-([a-z]+?){suffix_long}$", n, re.I
        )
        if m:
            locants = m.group(1)
            stem = m.group(2)
            return f"{stem}an-{locants}-{suffix_short}"

    # Alkenes / alkynes: "1-butene" → "but-1-ene", "1-pentyne" → "pent-1-yne"
    for unsuf in ["ene", "yne"]:
        m = re.match(rf"^(\d+(?:,\d+)*)-([a-z]+?){unsuf}$", n, re.I)
        if m:
            locants = m.group(1)
            stem = m.group(2)
            return f"{stem}-{locants}-{unsuf}"

    # Polyols: "1,2,3-propanetriol" → "propane-1,2,3-triol"
    m = re.match(r"^(\d+(?:,\d+)*)-([a-z]+?)(di|tri|tetra)(ol|amine|one)$", n, re.I)
    if m:
        locants, stem, multi, suf = m.groups()
        return f"{stem}-{locants}-{multi}{suf}"

    return name


def _ib_normalize(name: str) -> str:
    """Full IB normalization pipeline on a single candidate name."""
    result = _apply_trivial_mapping(name)
    result = _fix_old_locant_suffix(result)
    return result


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  LAYER 4 — IB-aware scoring for candidates                             ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def _ib_score(name: str) -> int:
    """
    Score a name for how well it conforms to IB Chemistry naming conventions.
    Higher is better.
    """
    n = name.lower().strip()
    score = 0

    # ── Strong IB substitutive patterns ──
    # Carboxylic acids: methanoic acid, ethanoic acid, propanoic acid …
    if re.search(r"(methan|ethan|propan|butan|pentan|hexan|heptan|octan)oic acid$", n):
        score += 10000
    # Alcohols: methanol, ethanol, propan-1-ol, butan-2-ol …
    if re.search(r"[a-z]+an-\d+(,\d+)*-ol$", n):
        score += 9800
    if re.match(r"^(methanol|ethanol)$", n):
        score += 9900
    # Aldehydes: methanal, ethanal, propanal …
    if re.match(r"^(methanal|ethanal|propanal|butanal|pentanal|hexanal)$", n):
        score += 9700
    # Ketones: propanone, butanone, pentan-2-one …
    if re.match(r"^(propanone|butanone)$", n):
        score += 9700
    if re.search(r"[a-z]+an-\d+-one$", n):
        score += 9600
    # Amines: methanamine, ethanamine, propan-1-amine …
    if re.search(r"[a-z]+an-\d+-amine$", n) or re.match(
        r"^(methanamine|ethanamine)$", n
    ):
        score += 9500
    # Amides: methanamide, ethanamide …
    if re.match(r"^[a-z]+(an)?amide$", n):
        score += 9400
    # Esters: methyl ethanoate, ethyl methanoate …
    if re.search(r"(methyl|ethyl|propyl|butyl) [a-z]+anoate$", n):
        score += 9300
    # Alkenes with modern locant: but-1-ene, pent-2-ene …
    if re.search(r"[a-z]+-\d+-ene$", n) or re.match(r"^(ethene|propene)$", n):
        score += 9200
    # Alkynes with modern locant: but-1-yne, pent-2-yne …
    if re.search(r"[a-z]+-\d+-yne$", n) or re.match(r"^(ethyne|propyne)$", n):
        score += 9100
    # Simple alkanes
    if re.match(
        r"^(methane|ethane|propane|butane|pentane|hexane|heptane|octane|nonane|decane)$",
        n,
    ):
        score += 9000
    # Cycloalkanes
    if re.match(r"^cyclo(propane|butane|pentane|hexane)$", n):
        score += 8900
    # IB-accepted aromatic retained names
    if n in ("benzene", "phenol", "benzoic acid", "phenylamine", "methylbenzene"):
        score += 9500
    # Branched alkane systematic: 2-methylpropane, 2,2-dimethylbutane …
    if re.match(r"^\d+(,\d+)*-(di|tri|tetra)?(methyl|ethyl|propyl)", n):
        score += 8500
    # Halogenoalkanes: 1-chloropropane, 2-bromobutane, chloromethane …
    if re.search(r"(chloro|bromo|iodo|fluoro)[a-z]+(ane|ene|yne)$", n):
        score += 8400
    if re.match(r"^\d+-?(chloro|bromo|iodo|fluoro)", n):
        score += 8300
    # Amino acids IB-style: 2-aminoethanoic acid
    if re.search(r"\d+-amino[a-z]+oic acid$", n):
        score += 9600
    # Diols/triols: propane-1,2,3-triol, butane-1,4-diol
    if re.search(r"[a-z]+ane-[\d,]+-[a-z]*ol$", n):
        score += 8800
    # Diacids: ethanedioic acid, propanedioic acid
    if re.search(r"[a-z]+anedioic acid$", n):
        score += 9200

    # ── Penalties ──
    # Trivial retained names that IB does NOT use
    trivial_bad = {
        "acetic acid", "formic acid", "propionic acid", "butyric acid",
        "formaldehyde", "acetaldehyde", "acetone", "ethylene", "propylene",
        "acetylene", "toluene", "aniline", "glycerol", "glycerine", "glycine",
        "alanine", "oxalic acid", "malonic acid", "succinic acid",
    }
    if n in trivial_bad:
        score -= 5000
    # Old-style locant (1-butanol instead of butan-1-ol)
    if re.match(r"^\d+-[a-z]+(ol|one|ene|yne|amine)$", n):
        score -= 2000
    # "X alcohol" style
    if re.search(r"\b(methyl|ethyl|propyl|butyl|isopropyl|isobutyl|tert-butyl) alcohol\b", n):
        score -= 3000
    # Noise
    if "%" in n or "solution" in n or "plastic container" in n or "glacial" in n:
        score -= 8000
    if len(n) > 90:
        score -= 2000

    return score


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Network helpers                                                        ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def _canonical_smiles(smiles: str) -> Optional[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def _fetch_pubchem_synonyms(cid: int) -> list[str]:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": _HTTP_UA})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
        block = data.get("InformationList", {}).get("Information", [{}])[0]
        return list(block.get("Synonym", []) or [])
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError, KeyError):
        return []


def _fetch_cactus_names(smiles: str) -> list[str]:
    q = urllib.parse.quote(smiles, safe="")
    url = f"https://cactus.nci.nih.gov/chemical/structure/{q}/names"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": _HTTP_UA})
        with urllib.request.urlopen(req, timeout=30) as resp:
            text = resp.read().decode(errors="replace")
        return [ln.strip() for ln in text.splitlines() if ln.strip()]
    except (urllib.error.URLError, urllib.error.HTTPError):
        return []


_GUNK_SYNONYM = re.compile(
    r"^[\d\-\.\s]+$|inchi=|^inchikey=|^inchi\s*=|chebi:|^nciopen|^wlns?:|^[A-Z]{1,3}\d{5,}_|_FLUKA$|_SIAL$|_ALDRICH$",
    re.I,
)


def _is_plausible_name_candidate(s: str) -> bool:
    t = s.strip()
    if len(t) < 3 or len(t) > 130:
        return False
    if not re.search(r"[a-zA-Z]", t):
        return False
    if _GUNK_SYNONYM.search(t):
        return False
    return True


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Master name selection: Dictionary → Normalize → Score                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def _preferred_display_name(
    target_smi: str,
    lexichem_name: str,
    user_query: str,
) -> str:
    """
    Pick the best IB-compliant name.
    Priority: 1) Hard-coded IB table  2) Lexichem + IB corrections.
    No extra network calls — fast.
    """
    # LAYER 1: exact match in IB dictionary (covers all common IB compounds)
    if target_smi and target_smi in _IB_TABLE:
        return _IB_TABLE[target_smi]

    # LAYER 2: take Lexichem name, apply trivial→systematic + locant fix
    base = lexichem_name or user_query or ""
    if base:
        return _ib_normalize(base)

    return user_query or ""


def lookup_iupac(name: str) -> CompoundInfo:
    """
    Resolve a name (common, trade, or IUPAC) via PubChem.
    Display name favors IUPAC substitutive-style strings when present in synonym lists.
    """
    query = name.strip()
    results = pcp.get_compounds(query, "name")
    if not results:
        raise ValueError(f"未找到名为 '{query}' 的化合物，请检查名称是否正确。")

    compound = results[0]
    lexichem = (compound.iupac_name or "").strip()
    raw_smi = (
        getattr(compound, "connectivity_smiles", None)
        or getattr(compound, "canonical_smiles", None)
        or ""
    )
    target_smi = _canonical_smiles(raw_smi) or ""

    preferred = _preferred_display_name(target_smi, lexichem, query)

    info = CompoundInfo(
        iupac_name=preferred,
        original_query=query,
        pubchem_lexichem_name=lexichem,
        molecular_formula=compound.molecular_formula or "",
        molecular_weight=compound.molecular_weight or 0.0,
        canonical_smiles=raw_smi or "",
        isomeric_smiles=(
            getattr(compound, "smiles", None)
            or getattr(compound, "isomeric_smiles", None)
            or ""
        ),
        inchi=compound.inchi or "",
        inchi_key=compound.inchikey or "",
        exact_mass=compound.exact_mass or 0.0,
        xlogp=compound.xlogp,
        hbond_donor_count=compound.h_bond_donor_count or 0,
        hbond_acceptor_count=compound.h_bond_acceptor_count or 0,
        rotatable_bond_count=compound.rotatable_bond_count or 0,
        heavy_atom_count=compound.heavy_atom_count or 0,
        charge=compound.charge or 0,
        cid=compound.cid or 0,
    )

    if compound.synonyms:
        info.synonyms = compound.synonyms[:10]
        if compound.synonyms:
            info.common_name = compound.synonyms[0]

    return info


def lookup_smiles(smiles: str) -> CompoundInfo:
    """
    Resolve a SMILES string (e.g. from a drawing editor) to full CompoundInfo.
    1. Validate & canonicalize with RDKit.
    2. Look up in PubChem by SMILES.
    3. Apply same systematic-name ranking as lookup_iupac.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("输入的结构无法被解析为有效分子，请检查画板中的结构。")

    canon = Chem.MolToSmiles(mol, canonical=True)

    results = pcp.get_compounds(canon, "smiles")
    if not results:
        return _build_info_from_rdkit_only(mol, canon)

    compound = results[0]
    lexichem = (compound.iupac_name or "").strip()
    raw_smi = (
        getattr(compound, "connectivity_smiles", None)
        or getattr(compound, "canonical_smiles", None)
        or canon
    )
    target_smi = _canonical_smiles(raw_smi) or canon

    preferred = _preferred_display_name(target_smi, lexichem, "")

    info = CompoundInfo(
        iupac_name=preferred or lexichem or canon,
        original_query=smiles,
        pubchem_lexichem_name=lexichem,
        molecular_formula=compound.molecular_formula or rdMolDescriptors.CalcMolFormula(mol),
        molecular_weight=compound.molecular_weight or Descriptors.ExactMolWt(mol),
        canonical_smiles=raw_smi or canon,
        isomeric_smiles=(
            getattr(compound, "smiles", None)
            or getattr(compound, "isomeric_smiles", None)
            or ""
        ),
        inchi=compound.inchi or (Chem.MolToInchi(mol) or ""),
        inchi_key=compound.inchikey or "",
        exact_mass=compound.exact_mass or Descriptors.ExactMolWt(mol),
        xlogp=compound.xlogp,
        hbond_donor_count=compound.h_bond_donor_count or 0,
        hbond_acceptor_count=compound.h_bond_acceptor_count or 0,
        rotatable_bond_count=compound.rotatable_bond_count or 0,
        heavy_atom_count=compound.heavy_atom_count or mol.GetNumHeavyAtoms(),
        charge=compound.charge or Chem.GetFormalCharge(mol),
        cid=compound.cid or 0,
    )
    if compound.synonyms:
        info.synonyms = compound.synonyms[:10]
        if compound.synonyms:
            info.common_name = compound.synonyms[0]
    return info


def _build_info_from_rdkit_only(mol, canon: str) -> CompoundInfo:
    """Fallback when PubChem has no match: compute everything locally via RDKit."""
    inchi = Chem.MolToInchi(mol) or ""
    inchi_key = Chem.InchiToInchiKey(inchi) if inchi else ""
    formula = rdMolDescriptors.CalcMolFormula(mol)

    # IB table first
    if canon in _IB_TABLE:
        best = _IB_TABLE[canon]
    else:
        cactus_names = _fetch_cactus_names(canon)
        candidates = [n for n in cactus_names if _is_plausible_name_candidate(n)]
        if candidates:
            pairs = [(_ib_normalize(c), _ib_score(_ib_normalize(c))) for c in candidates]
            pairs.sort(key=lambda x: (-x[1], len(x[0])))
            best = pairs[0][0] if pairs[0][1] > 0 else canon
        else:
            best = canon

    return CompoundInfo(
        iupac_name=best,
        original_query=canon,
        molecular_formula=formula,
        molecular_weight=Descriptors.ExactMolWt(mol),
        canonical_smiles=canon,
        inchi=inchi,
        inchi_key=inchi_key or "",
        exact_mass=Descriptors.ExactMolWt(mol),
        heavy_atom_count=mol.GetNumHeavyAtoms(),
        charge=Chem.GetFormalCharge(mol),
    )


def validate_structure(smiles: str) -> tuple[bool, str]:
    """
    Basic chemical sanity check on a drawn structure.
    Returns (is_valid, message).
    """
    if not smiles or not smiles.strip():
        return False, "画板为空，请先绘制一个分子结构。"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "无法解析该结构，可能含有无效键或不合法的原子。"

    problems = Chem.DetectChemistryProblems(mol)
    if problems:
        msgs = "; ".join(p.Message() for p in problems)
        return False, f"结构存在化学问题: {msgs}"

    num_atoms = mol.GetNumHeavyAtoms()
    if num_atoms == 0:
        return False, "分子中没有重原子。"
    if num_atoms > 200:
        return False, "分子过大（>200 重原子），请简化。"

    return True, "结构有效 ✓"


def smiles_to_mol(smiles: str):
    """Convert SMILES to an RDKit Mol object."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法从SMILES '{smiles}' 生成分子对象。")
    return mol


def generate_2d_image(smiles: str, size: tuple = (500, 400)) -> bytes:
    """Generate a 2D structure PNG image from SMILES, returned as bytes."""
    mol = smiles_to_mol(smiles)
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def get_empirical_formula(smiles: str) -> str:
    """Return the molecular formula computed by RDKit."""
    mol = smiles_to_mol(smiles)
    return rdMolDescriptors.CalcMolFormula(mol)


def get_molecular_weight(smiles: str) -> float:
    mol = smiles_to_mol(smiles)
    return Descriptors.ExactMolWt(mol)


def smiles_to_inchi(smiles: str) -> str:
    mol = smiles_to_mol(smiles)
    return Chem.MolToInchi(mol) or ""


def smiles_to_inchi_key(smiles: str) -> str:
    inchi = smiles_to_inchi(smiles)
    if inchi:
        return Chem.InchiToInchiKey(inchi) or ""
    return ""


def format_molecular_formula_html(formula: str) -> str:
    """Convert a molecular formula like C2H6O into HTML with subscripts."""
    html = ""
    for ch in formula:
        if ch.isdigit():
            html += f"<sub>{ch}</sub>"
        else:
            html += ch
    return html


def generate_3d_mol_block(smiles: str) -> str:
    """Generate 3D-optimised mol block for interactive viewing."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    res = AllChem.EmbedMolecule(mol, params)
    if res == -1:
        params.useRandomCoords = True
        AllChem.EmbedMolecule(mol, params)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            pass
    return Chem.MolToMolBlock(mol)


def search_formula_isomers(formula: str, max_results: int = 25) -> list[dict]:
    """
    Search PubChem for structural isomers sharing *formula*.
    Returns list of dicts with CID, CanonicalSMILES, IUPACName,
    MolecularWeight, and iupac_ib (IB-normalised display name).
    """
    formula = formula.strip()
    if not formula or not re.match(r"^[A-Z][A-Za-z0-9]*$", formula):
        raise ValueError("请输入有效的分子式，如 C4H10O")

    encoded = urllib.parse.quote(formula, safe="")
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    url = f"{base}/compound/formula/{encoded}/cids/JSON"
    req = urllib.request.Request(url, headers={"User-Agent": _HTTP_UA})

    def _pug_get(u, retries=3):
        """GET with retry to survive transient SSL / connection drops."""
        for attempt in range(retries):
            r = urllib.request.Request(u, headers={"User-Agent": _HTTP_UA})
            try:
                with urllib.request.urlopen(r, timeout=30) as resp:
                    return json.loads(resp.read().decode())
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    return None
                if attempt == retries - 1:
                    raise ValueError(f"PubChem error: HTTP {e.code}")
            except urllib.error.URLError:
                if attempt == retries - 1:
                    return None
            time.sleep(1)
        return None

    data = _pug_get(url)
    if data is None:
        return []

    retries = 0
    while "Waiting" in data and retries < 20:
        list_key = data["Waiting"]["ListKey"]
        time.sleep(2)
        poll_url = f"{base}/compound/listkey/{list_key}/cids/JSON"
        result = _pug_get(poll_url)
        if result is not None:
            data = result
        retries += 1

    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return []

    batch = cids[: min(len(cids), 100)]
    cid_str = ",".join(str(c) for c in batch)
    prop_url = (
        f"{base}/compound/cid/{cid_str}/property/"
        "IUPACName,CanonicalSMILES,ConnectivitySMILES,MolecularFormula,MolecularWeight/JSON"
    )
    prop_data = _pug_get(prop_url)
    if prop_data is None:
        return []

    props = prop_data.get("PropertyTable", {}).get("Properties", [])

    seen: set[str] = set()
    results: list[dict] = []
    for p in props:
        smi = p.get("CanonicalSMILES") or p.get("ConnectivitySMILES") or ""
        if not smi or smi in seen:
            continue
        seen.add(smi)
        canon = _canonical_smiles(smi) or smi
        lexichem = p.get("IUPACName", "")
        ib_name = _preferred_display_name(canon, lexichem, "")
        p["iupac_ib"] = ib_name or lexichem or smi
        results.append(p)
        if len(results) >= max_results:
            break

    return results


def get_condensed_formula(smiles: str) -> str:
    """
    Attempt to build a condensed structural formula from SMILES.
    For simple chains this works well (e.g. CH3CH2OH).
    Falls back to canonical SMILES for complex molecules.
    """
    mol = smiles_to_mol(smiles)
    mol = Chem.AddHs(mol)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return Chem.MolToSmiles(Chem.RemoveHs(mol))

    atoms = mol.GetAtoms()
    if mol.GetNumAtoms() == 0:
        return smiles

    visited = set()
    parts = []

    start_idx = None
    for atom in atoms:
        if atom.GetDegree() == 1 and atom.GetSymbol() != "H":
            start_idx = atom.GetIdx()
            break
    if start_idx is None:
        for atom in atoms:
            if atom.GetSymbol() != "H":
                start_idx = atom.GetIdx()
                break
    if start_idx is None:
        return smiles

    def _build(idx):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == "H":
            return
        visited.add(idx)

        h_count = sum(
            1
            for n in atom.GetNeighbors()
            if n.GetSymbol() == "H" and n.GetIdx() not in visited
        )
        for n in atom.GetNeighbors():
            if n.GetSymbol() == "H":
                visited.add(n.GetIdx())

        symbol = atom.GetSymbol()
        if h_count > 1:
            symbol += f"H{h_count}"
        elif h_count == 1:
            symbol += "H"

        parts.append(symbol)

        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited and neighbor.GetSymbol() != "H":
                bond = mol.GetBondBetweenAtoms(idx, neighbor.GetIdx())
                if bond.GetBondTypeAsDouble() == 2.0:
                    parts.append("=")
                elif bond.GetBondTypeAsDouble() == 3.0:
                    parts.append("≡")
                _build(neighbor.GetIdx())

    _build(start_idx)
    result = "".join(parts)
    return result if result else smiles
