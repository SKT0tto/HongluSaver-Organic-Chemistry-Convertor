"""
HongluSaver - Chemistry formula converter
Author: Henry Lin
"""

from dataclasses import replace
import json as _json

import streamlit as st
import streamlit.components.v1 as components
from streamlit_ketcher import st_ketcher
from chemistry_engine import (
    lookup_iupac,
    lookup_smiles,
    validate_structure,
    generate_2d_image,
    generate_3d_mol_block,
    search_formula_isomers,
    get_condensed_formula,
    format_molecular_formula_html,
)

# ── i18n ─────────────────────────────────────────────────────────────────────

_T = {
    "zh": {
        "page_title_cn": "🧪 IUPAC 有机化学命名 → 结构转换器",
        "page_desc": "名称查询 · 结构绘制 · 多种化学表示一键生成",
        "quick_examples": "📋 快速示例（查询模式）",
        "click_to_fill": "点击按钮将名称填入输入框",
        "cache_note": "重复查询自动走本地缓存（24 h 过期）。",
        "data_source": "数据来源",
        "structure_render": "结构渲染",
        "editor": "画板编辑器",
        "tab_search": "🔍 查询模式",
        "tab_draw": "🎨 自由绘制模式",
        "input_hint": "##### 输入化合物名称（俗名 / 商品名 / IUPAC 均可）",
        "input_placeholder": "例如: acetic acid、ethanol、2-methylpropane、aspirin …",
        "btn_search": "🔍 查询",
        "searching": "正在查询 **{q}** …",
        "warn_empty": "⚠️ 请输入一个有效的化合物名称（英文）。",
        "draw_title": "##### 🖊️ 在下方画板绘制 Skeletal Formula，然后点击「识别结构」",
        "draw_hint": "使用 Ketcher 编辑器：左侧工具栏选择原子（C, N, O …）和键类型，直接在画布上拖拽绘制。完成后点下方按钮。",
        "btn_analyze": "🧬 识别结构 → 生成 IUPAC 命名",
        "warn_empty_canvas": "画板为空，请先绘制一个分子结构再点击识别。",
        "analyzing": "正在解析并查找 IUPAC 命名…",
        "parse_fail": "❌ 解析失败: {e}",
        "you_entered": "您输入: **{q}**",
        "lexichem_label": "PubChem Lexichem 名",
        "common_name_label": "常用名 / 同义词",
        "struct_2d": "🔬 2D 结构图",
        "struct_fail": "无法生成结构图: {e}",
        "chem_repr": "📊 化学表示",
        "mol_formula": "分子式 Molecular Formula",
        "condensed": "简缩结构式 Condensed Structural Formula",
        "isomeric_smiles": "异构 SMILES (Isomeric)",
        "mol_props": "⚗️ 分子属性",
        "mw": "分子量",
        "exact_mass": "精确质量",
        "heavy_atoms": "重原子数",
        "formal_charge": "形式电荷",
        "hbd": "氢键供体数",
        "hba": "氢键受体数",
        "rotatable": "可旋转键数",
        "synonyms_title": "🏷️ 其他名称 / 同义词",
        "viewer_3d": "🧊 3D 分子模型",
        "viewer_3d_hint": "拖拽旋转 · 滚轮缩放 · 自动旋转中",
        "3d_fail": "无法生成 3D 模型: {e}",
        "tab_isomer": "🧩 同分异构体",
        "isomer_input_hint": "##### 输入分子式，查找所有结构异构体",
        "isomer_placeholder": "例如: C4H10O, C3H8O, C2H6O …",
        "btn_find_isomers": "🔍 查找",
        "isomer_searching": "正在搜索 **{f}** 的同分异构体…",
        "isomer_found": "找到 **{n}** 个结构不同的化合物",
        "isomer_none": "⚠️ 未找到该分子式对应的化合物，请检查分子式。",
        "isomer_invalid": "⚠️ 请输入有效的分子式（如 C4H10O）。",
        "isomer_examples": "常见 IB 分子式",
        "footer": "HongluSaver — Chemistry formula converter · 作者 Author: Henry Lin",
    },
    "en": {
        "page_title_cn": "🧪 IUPAC Organic Name → Structure Converter",
        "page_desc": "Name lookup · Structure drawing · Multiple chemical representations",
        "quick_examples": "📋 Quick Examples (Search)",
        "click_to_fill": "Click a button to fill the input",
        "cache_note": "Repeated queries use local cache (expires 24 h).",
        "data_source": "Data source",
        "structure_render": "Structure rendering",
        "editor": "Drawing editor",
        "tab_search": "🔍 Search Mode",
        "tab_draw": "🎨 Free Drawing Mode",
        "input_hint": "##### Enter compound name (common / trade / IUPAC)",
        "input_placeholder": "e.g. acetic acid, ethanol, 2-methylpropane, aspirin …",
        "btn_search": "🔍 Search",
        "searching": "Searching **{q}** …",
        "warn_empty": "⚠️ Please enter a valid compound name.",
        "draw_title": "##### 🖊️ Draw a skeletal formula below, then click 'Identify'",
        "draw_hint": "Use the Ketcher editor: select atoms (C, N, O …) and bond types from the left toolbar, drag on the canvas to draw. Click the button below when done.",
        "btn_analyze": "🧬 Identify Structure → Generate IUPAC Name",
        "warn_empty_canvas": "Canvas is empty. Please draw a molecular structure first.",
        "analyzing": "Parsing and looking up IUPAC name…",
        "parse_fail": "❌ Parse failed: {e}",
        "you_entered": "You entered: **{q}**",
        "lexichem_label": "PubChem Lexichem name",
        "common_name_label": "Common name / Synonym",
        "struct_2d": "🔬 2D Structure",
        "struct_fail": "Cannot generate structure image: {e}",
        "chem_repr": "📊 Chemical Representations",
        "mol_formula": "Molecular Formula",
        "condensed": "Condensed Structural Formula",
        "isomeric_smiles": "Isomeric SMILES",
        "mol_props": "⚗️ Molecular Properties",
        "mw": "Molecular Weight",
        "exact_mass": "Exact Mass",
        "heavy_atoms": "Heavy Atoms",
        "formal_charge": "Formal Charge",
        "hbd": "H-Bond Donors",
        "hba": "H-Bond Acceptors",
        "rotatable": "Rotatable Bonds",
        "synonyms_title": "🏷️ Other Names / Synonyms",
        "viewer_3d": "🧊 3D Model",
        "viewer_3d_hint": "Drag to rotate · Scroll to zoom · Auto-spinning",
        "3d_fail": "Cannot generate 3D model: {e}",
        "tab_isomer": "🧩 Isomers",
        "isomer_input_hint": "##### Enter a molecular formula to find structural isomers",
        "isomer_placeholder": "e.g. C4H10O, C3H8O, C2H6O …",
        "btn_find_isomers": "🔍 Search",
        "isomer_searching": "Searching structural isomers of **{f}**…",
        "isomer_found": "Found **{n}** structurally different compounds",
        "isomer_none": "⚠️ No compounds found for this formula. Check the formula.",
        "isomer_invalid": "⚠️ Please enter a valid molecular formula (e.g. C4H10O).",
        "isomer_examples": "Common IB formulas",
        "footer": "HongluSaver — Chemistry formula converter · Author: Henry Lin",
    },
}


def t(key: str, **kwargs) -> str:
    lang = st.session_state.get("lang", "zh")
    text = _T.get(lang, _T["zh"]).get(key, key)
    if kwargs:
        text = text.format(**kwargs)
    return text


# ── Cached helpers ───────────────────────────────────────────────────────────

@st.cache_data(ttl=86_400, max_entries=320, show_spinner=False)
def _cached_compound_lookup(query_normalized: str):
    return lookup_iupac(query_normalized)


@st.cache_data(ttl=86_400, max_entries=180, show_spinner=False)
def _cached_structure_png(canonical_smiles: str, width: int, height: int) -> bytes:
    return generate_2d_image(canonical_smiles, size=(width, height))


@st.cache_data(ttl=86_400, max_entries=320, show_spinner=False)
def _cached_smiles_lookup(canonical_smiles: str):
    return lookup_smiles(canonical_smiles)


@st.cache_data(ttl=86_400, max_entries=200, show_spinner=False)
def _cached_3d_mol_block(canonical_smiles: str) -> str:
    return generate_3d_mol_block(canonical_smiles)


@st.cache_data(ttl=86_400, max_entries=100, show_spinner=False)
def _cached_isomer_search(formula: str):
    return search_formula_isomers(formula)


# ── Page Config ──────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="HongluSaver - Chemistry formula converter",
    page_icon="🧪",
    layout="wide",
)

# ── Custom CSS ───────────────────────────────────────────────────────────────

st.markdown(
    """
<style>
    @import url('https://fonts.googleapis.com/css2?family=Exo+2:wght@500;600;700&family=Noto+Sans+SC:wght@600;700&family=Orbitron:wght@600;700;800;900&display=swap');

    @keyframes brand-shine {
        0%, 100% { background-position: 0% 50%; }
        50% { background-position: 100% 50%; }
    }

    .block-container {
        padding-top: 2rem !important;
        overflow-x: visible !important;
    }
    section[data-testid="stMain"] > div {
        overflow: visible !important;
    }
    .main-header {
        text-align: center;
        padding: 1.5rem 0.5rem 1rem 0.5rem;
        overflow: visible;
    }
    .main-header .brand-shell {
        display: inline-block;
        max-width: 100%;
        padding: 0.25em 0.15em 0.45em;
        margin: 0 auto;
        box-sizing: border-box;
        overflow: visible;
    }
    .main-header .brand-name {
        font-family: 'Orbitron', 'Segoe UI', system-ui, sans-serif;
        font-size: clamp(2.75rem, 7.5vw, 4.35rem);
        font-weight: 900;
        letter-spacing: 0.06em;
        margin: 0;
        line-height: 1.22;
        padding: 0.08em 0;
        background: linear-gradient(105deg, #7c3aed 0%, #6366f1 35%, #a855f7 55%, #ec4899 100%);
        background-size: 200% auto;
        -webkit-background-clip: text;
        background-clip: text;
        -webkit-text-fill-color: transparent;
        animation: brand-shine 10s ease-in-out infinite;
    }
    .main-header .brand-tagline {
        font-family: 'Exo 2', 'Noto Sans SC', sans-serif;
        color: #4b5563;
        font-size: clamp(1.05rem, 2.2vw, 1.45rem);
        font-weight: 600;
        letter-spacing: 0.18em;
        text-transform: uppercase;
        margin: 0.35rem 0 0.45rem 0;
    }
    .main-header .author-credit {
        color: #999;
        font-size: 0.9rem;
        margin: 0 0 0.75rem 0;
    }
    .main-header .page-title-cn {
        font-family: 'Noto Sans SC', 'PingFang SC', 'Microsoft YaHei', sans-serif;
        color: #333;
        font-size: clamp(1.2rem, 2.8vw, 1.65rem);
        font-weight: 700;
        margin: 0.5rem 0 0.35rem 0;
    }
    .main-header p.desc {
        font-family: 'Noto Sans SC', 'PingFang SC', 'Microsoft YaHei', sans-serif;
        color: #888;
        font-size: 1.05rem;
        margin: 0;
    }
    .formula-card {
        background: #f8f9fc;
        border-radius: 12px;
        padding: 1.2rem 1.4rem;
        margin-bottom: 0.8rem;
        border: 1px solid #e4e7f0;
    }
    .formula-card h4 {
        margin: 0 0 0.3rem 0;
        color: #555;
        font-size: 0.85rem;
        text-transform: uppercase;
        letter-spacing: 0.04em;
    }
    .formula-card .value {
        font-size: 1.2rem;
        font-weight: 600;
        color: #1a1a2e;
        word-break: break-all;
    }
    .struct-img {
        display: flex;
        justify-content: center;
        padding: 1rem;
    }
    div[data-testid="stMetric"] {
        background: #f8f9fc;
        border: 1px solid #e4e7f0;
        border-radius: 12px;
        padding: 1rem;
    }

    /* ── Sidebar brand logo ── */
    .sidebar-logo {
        text-align: center;
        padding: 0.6rem 0 0.3rem 0;
    }
    .sidebar-logo .logo-benzene {
        display: inline-block;
        position: relative;
        width: 120px;
        height: 120px;
        margin-bottom: 0.35rem;
    }
    .sidebar-logo .logo-benzene svg {
        width: 100%;
        height: 100%;
        filter: drop-shadow(0 3px 10px rgba(124, 58, 237, 0.25));
    }
    @media (prefers-color-scheme: dark) {
        .sidebar-logo .logo-benzene .logo-txt-stroke {
            stroke: #1a1a2e !important;
        }
    }
    [data-theme="dark"] .sidebar-logo .logo-benzene .logo-txt-stroke,
    .stApp[data-theme="dark"] .sidebar-logo .logo-benzene .logo-txt-stroke {
        stroke: #1a1a2e !important;
    }
    .sidebar-logo .logo-sub {
        font-family: 'Exo 2', sans-serif;
        font-size: 0.72rem;
        font-weight: 600;
        letter-spacing: 0.14em;
        text-transform: uppercase;
        color: #6b7280;
        margin: 0.1rem 0 0 0;
    }
    .sidebar-logo .logo-author {
        font-size: 0.7rem;
        color: #aaa;
        margin: 0.35rem 0 0 0;
    }
</style>
""",
    unsafe_allow_html=True,
)

# ── Header ───────────────────────────────────────────────────────────────────

st.markdown(
    f"""
<div class="main-header">
    <div class="brand-shell">
    <p class="brand-name">HongluSaver</p>
    </div>
    <p class="brand-tagline">Chemistry formula converter</p>
    <p class="author-credit">Made by Henry Lin</p>
    <p class="page-title-cn">{t("page_title_cn")}</p>
    <p class="desc">{t("page_desc")}</p>
</div>
""",
    unsafe_allow_html=True,
)

# ── Sidebar ──────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown(
        """
<div class="sidebar-logo">
    <div class="logo-benzene">
        <svg viewBox="0 0 120 120" xmlns="http://www.w3.org/2000/svg">
            <defs>
                <linearGradient id="grad1" x1="0%" y1="0%" x2="100%" y2="100%">
                    <stop offset="0%"   stop-color="#7c3aed"/>
                    <stop offset="50%"  stop-color="#a855f7"/>
                    <stop offset="100%" stop-color="#ec4899"/>
                </linearGradient>
            </defs>
            <!-- outer hexagon -->
            <polygon points="60,10 103.3,35 103.3,85 60,110 16.7,85 16.7,35"
                     fill="none" stroke="url(#grad1)" stroke-width="2.8"
                     stroke-linejoin="round"/>
            <!-- inner circle (aromatic ring) -->
            <circle cx="60" cy="60" r="33" fill="none"
                    stroke="url(#grad1)" stroke-width="2.8"/>
            <!-- brand text: stroke layer (outline) -->
            <text class="logo-txt-stroke" x="60" y="62"
                  text-anchor="middle" dominant-baseline="central"
                  transform="rotate(-18 60 60)"
                  font-family="'Orbitron','Segoe UI',system-ui,sans-serif"
                  font-size="13" font-weight="900" letter-spacing="0.4"
                  fill="none" stroke="#ffffff" stroke-width="4"
                  stroke-linejoin="round">HongluSaver</text>
            <!-- brand text: fill layer (gradient on top) -->
            <text x="60" y="62"
                  text-anchor="middle" dominant-baseline="central"
                  transform="rotate(-18 60 60)"
                  font-family="'Orbitron','Segoe UI',system-ui,sans-serif"
                  font-size="13" font-weight="900" letter-spacing="0.4"
                  fill="url(#grad1)">HongluSaver</text>
        </svg>
    </div>
    <p class="logo-sub">Chemistry Formula Converter</p>
    <p class="logo-author">by Henry Lin</p>
</div>
""",
        unsafe_allow_html=True,
    )

    # ── Language toggle ──
    lang_options = {"zh": "🇨🇳 中文", "en": "🇬🇧 English"}
    current = st.session_state.get("lang", "zh")
    chosen = st.radio(
        "Language",
        options=list(lang_options.keys()),
        format_func=lambda k: lang_options[k],
        index=list(lang_options.keys()).index(current),
        key="lang_radio",
        horizontal=True,
        label_visibility="collapsed",
    )
    if chosen != current:
        st.session_state["lang"] = chosen
        st.rerun()

    st.divider()
    st.header(t("quick_examples"))
    st.caption(t("click_to_fill"))
    examples = [
        "ethanol",
        "acetic acid",
        "benzene",
        "acetone",
        "propan-1-ol",
        "butanoic acid",
        "cyclohexane",
        "2-methylpropane",
        "phenol",
        "glycine",
        "aspirin",
        "caffeine",
        "glucose",
        "toluene",
        "aniline",
    ]
    for ex in examples:
        if st.button(ex, key=f"ex_{ex}", use_container_width=True):
            st.session_state["search_name_input"] = ex
            st.session_state["auto_search"] = True
            st.rerun()

    st.divider()
    st.caption(t("cache_note"))
    st.markdown(
        f"**{t('data_source')}**: [PubChem](https://pubchem.ncbi.nlm.nih.gov/)  \n"
        f"**{t('structure_render')}**: [RDKit](https://www.rdkit.org/)  \n"
        f"**{t('editor')}**: [Ketcher](https://lifescience.opensource.epam.com/ketcher/)"
    )


# ══════════════════════════════════════════════════════════════════════════════
# Helper: display a CompoundInfo block (shared between both tabs)
# ══════════════════════════════════════════════════════════════════════════════

def _render_3d_viewer(smiles: str, height: int = 420):
    """Render an interactive 3D molecule viewer using 3Dmol.js from CDN."""
    try:
        mol_block = _cached_3d_mol_block(smiles)
    except Exception as e:
        st.warning(t("3d_fail", e=str(e)))
        return

    mol_block_js = _json.dumps(mol_block)
    viewer_id = f"v{hash(smiles) & 0xFFFFFFFF}"

    html = f"""
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <div id="{viewer_id}"
         style="width:100%;height:{height}px;position:relative;
                border:1px solid #e4e7f0;border-radius:12px;overflow:hidden;">
    </div>
    <script>
    (function(){{
        var v=$3Dmol.createViewer("{viewer_id}",{{backgroundColor:"white"}});
        v.addModel({mol_block_js},"sdf");
        v.setStyle({{}},{{stick:{{radius:0.14,colorscheme:"Jmol"}},
                         sphere:{{scale:0.28,colorscheme:"Jmol"}}}});
        v.zoomTo();v.render();v.spin("y",0.8);
    }})();
    </script>
    """
    components.html(html, height=height + 10)


def _render_compound_result(info, *, show_query_caption: bool = True):
    st.markdown(f"## 📄 {info.iupac_name}")
    if show_query_caption and info.original_query and info.original_query.lower() != info.iupac_name.lower():
        st.caption(t("you_entered", q=info.original_query))
    if info.pubchem_lexichem_name and info.pubchem_lexichem_name.lower() != info.iupac_name.lower():
        st.caption(f"{t('lexichem_label')}: {info.pubchem_lexichem_name}")
    if info.common_name and info.common_name.lower() != info.iupac_name.lower():
        st.caption(f"{t('common_name_label')}: {info.common_name}")
    if info.cid:
        st.caption(f"PubChem CID: [{info.cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{info.cid})")

    st.divider()

    col_img, col_data = st.columns([1, 1], gap="large")

    with col_img:
        tab_2d, tab_3d = st.tabs([t("struct_2d"), t("viewer_3d")])
        with tab_2d:
            try:
                img_bytes = _cached_structure_png(info.canonical_smiles, 550, 450)
                st.image(img_bytes, use_container_width=True)
            except Exception as e:
                st.warning(t("struct_fail", e=e))
        with tab_3d:
            st.caption(t("viewer_3d_hint"))
            _render_3d_viewer(info.canonical_smiles)

    with col_data:
        st.subheader(t("chem_repr"))

        formula_html = format_molecular_formula_html(info.molecular_formula)
        st.markdown(
            f'<div class="formula-card"><h4>{t("mol_formula")}</h4>'
            f'<div class="value">{formula_html}</div></div>',
            unsafe_allow_html=True,
        )

        try:
            condensed = get_condensed_formula(info.canonical_smiles)
        except Exception:
            condensed = info.canonical_smiles
        st.markdown(
            f'<div class="formula-card"><h4>{t("condensed")}</h4>'
            f'<div class="value">{condensed}</div></div>',
            unsafe_allow_html=True,
        )

        st.markdown(
            f'<div class="formula-card"><h4>SMILES</h4>'
            f'<div class="value">{info.canonical_smiles}</div></div>',
            unsafe_allow_html=True,
        )

        if info.isomeric_smiles and info.isomeric_smiles != info.canonical_smiles:
            st.markdown(
                f'<div class="formula-card"><h4>{t("isomeric_smiles")}</h4>'
                f'<div class="value">{info.isomeric_smiles}</div></div>',
                unsafe_allow_html=True,
            )

        st.markdown(
            f'<div class="formula-card"><h4>InChI</h4>'
            f'<div class="value" style="font-size:0.95rem">{info.inchi}</div></div>',
            unsafe_allow_html=True,
        )

        st.markdown(
            f'<div class="formula-card"><h4>InChI Key</h4>'
            f'<div class="value">{info.inchi_key}</div></div>',
            unsafe_allow_html=True,
        )

    st.divider()

    st.subheader(t("mol_props"))
    m1, m2, m3, m4 = st.columns(4)
    m1.metric(t("mw"), f"{info.molecular_weight:.3f} g/mol")
    m2.metric(t("exact_mass"), f"{info.exact_mass:.4f}")
    m3.metric(t("heavy_atoms"), info.heavy_atom_count)
    m4.metric(t("formal_charge"), info.charge)

    m5, m6, m7, m8 = st.columns(4)
    m5.metric(t("hbd"), info.hbond_donor_count)
    m6.metric(t("hba"), info.hbond_acceptor_count)
    m7.metric(t("rotatable"), info.rotatable_bond_count)
    m8.metric("XLogP", f"{info.xlogp:.2f}" if info.xlogp is not None else "N/A")

    if info.synonyms:
        st.divider()
        st.subheader(t("synonyms_title"))
        st.write(" · ".join(info.synonyms))


# ══════════════════════════════════════════════════════════════════════════════
# Two-tab layout
# ══════════════════════════════════════════════════════════════════════════════

tab_search, tab_draw, tab_isomer = st.tabs([t("tab_search"), t("tab_draw"), t("tab_isomer")])

# ── Tab 1: Search by name ────────────────────────────────────────────────────

with tab_search:
    st.markdown(t("input_hint"))
    col_input, col_btn = st.columns([5, 1])
    with col_input:
        name = st.text_input(
            "compound name",
            placeholder=t("input_placeholder"),
            label_visibility="collapsed",
            key="search_name_input",
        )
    with col_btn:
        search = st.button(t("btn_search"), type="primary", use_container_width=True)

    auto = st.session_state.pop("auto_search", False)

    if (search or auto) and name.strip():
        query = name.strip()
        with st.spinner(t("searching", q=query)):
            try:
                base = _cached_compound_lookup(query.lower())
                info = replace(base, original_query=query)
            except Exception as e:
                st.error(f"❌ {e}")
                st.stop()
        _render_compound_result(info)
    elif search:
        st.warning(t("warn_empty"))

# ── Tab 2: Free drawing mode ─────────────────────────────────────────────────

with tab_draw:
    st.markdown(t("draw_title"))
    st.caption(t("draw_hint"))

    drawn_smiles = st_ketcher(value="", height=520, key="ketcher_editor")

    st.divider()

    analyze = st.button(t("btn_analyze"), type="primary", use_container_width=True)

    if analyze:
        if not drawn_smiles or not drawn_smiles.strip():
            st.warning(t("warn_empty_canvas"))
        else:
            ok, msg = validate_structure(drawn_smiles)
            if not ok:
                st.error(f"❌ {msg}")
            else:
                st.success(msg)
                with st.spinner(t("analyzing")):
                    try:
                        from rdkit import Chem
                        canon = Chem.MolToSmiles(
                            Chem.MolFromSmiles(drawn_smiles), canonical=True
                        )
                        info = _cached_smiles_lookup(canon)
                    except Exception as e:
                        st.error(t("parse_fail", e=e))
                        st.stop()
                _render_compound_result(info, show_query_caption=False)


# ── Tab 3: Isomer generator ───────────────────────────────────────────────

with tab_isomer:
    st.markdown(t("isomer_input_hint"))
    col_f, col_b = st.columns([5, 1])
    with col_f:
        formula_input = st.text_input(
            "formula",
            placeholder=t("isomer_placeholder"),
            label_visibility="collapsed",
            key="isomer_formula_input",
        )
    with col_b:
        search_iso = st.button(
            t("btn_find_isomers"), type="primary", use_container_width=True,
            key="btn_isomer_search",
        )

    def _set_iso_example(val):
        st.session_state["isomer_formula_input"] = val
        st.session_state["auto_iso_search"] = True

    st.caption(t("isomer_examples"))
    iso_examples = ["C4H10O", "C3H8O", "C2H6O", "C4H8", "C3H6O", "C4H10", "C3H6O2"]
    iso_cols = st.columns(len(iso_examples))
    for i, ef in enumerate(iso_examples):
        with iso_cols[i]:
            st.button(ef, key=f"iso_ex_{ef}", use_container_width=True,
                      on_click=_set_iso_example, args=(ef,))

    auto_iso = st.session_state.pop("auto_iso_search", False)

    if (search_iso or auto_iso) and formula_input.strip():
        with st.spinner(t("isomer_searching", f=formula_input.strip())):
            try:
                iso_results = _cached_isomer_search(formula_input.strip())
            except Exception as e:
                st.error(f"❌ {e}")
                st.stop()

        if not iso_results:
            st.warning(t("isomer_none"))
        else:
            st.success(t("isomer_found", n=len(iso_results)))
            for row_start in range(0, len(iso_results), 3):
                grid = st.columns(3)
                for j, col in enumerate(grid):
                    idx = row_start + j
                    if idx >= len(iso_results):
                        break
                    iso = iso_results[idx]
                    with col:
                        smi = iso.get("CanonicalSMILES") or iso.get("ConnectivitySMILES") or ""
                        name = iso.get("iupac_ib", smi)
                        try:
                            img = _cached_structure_png(smi, 300, 250)
                            st.image(img, use_container_width=True)
                        except Exception:
                            st.info("(structure unavailable)")
                        st.markdown(f"**{name}**")
                        mw = iso.get("MolecularWeight")
                        if mw:
                            st.caption(f"MW: {float(mw):.2f} · `{smi}`")
                        else:
                            st.caption(f"`{smi}`")
    elif search_iso:
        st.warning(t("isomer_invalid"))


# ── Footer ───────────────────────────────────────────────────────────────────

st.divider()
st.caption(t("footer"))
