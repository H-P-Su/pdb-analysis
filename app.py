#!/usr/bin/env python3
"""
app.py — PDB Toolkit multi-page Streamlit app.

Run:
    streamlit run app.py
"""

import streamlit as st

st.set_page_config(
    page_title="PDB Toolkit",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

pg = st.navigation([
    st.Page("app_structure.py",  title="Structure Analyzer", icon="🔬"),
    st.Page("app_alphafold.py",  title="AlphaFold + Variants", icon="🧬"),
])
pg.run()
