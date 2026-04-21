#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
差异代谢物药用筛选 Web 应用
===========================
基于 Streamlit 构建，支持文件上传、可视化、明星分子精选、KEGG 通路富集

用法：
    pip install -r requirements.txt
    streamlit run app.py
"""

import io
import os
import re
import time
import sys
from datetime import datetime

import pandas as pd
import numpy as np

# 新增模块
try:
    from report_generator import generate_word_report
except ImportError:
    generate_word_report = None

try:
    from pharma_db import (
        load_tcmsp, load_tcmsp_from_dir, load_drugbank, load_ttd,
        match_pharma_db, compute_pharma_evidence_score
    )
except ImportError:
    load_tcmsp = load_tcmsp_from_dir = load_drugbank = load_ttd = match_pharma_db = compute_pharma_evidence_score = None

try:
    from pharma_cache import match_pharma_online, query_pharma_info, compute_pharma_evidence_score as compute_pharma_evidence_score_online
except ImportError:
    match_pharma_online = query_pharma_info = compute_pharma_evidence_score_online = None

try:
    from multivariate_stats import (
        detect_species_columns_from_diff, build_abundance_matrix, compute_log2fc,
        run_pca, run_plsda, run_oplsda, permutation_test_plsda,
        plot_scores_scatter, plot_loading_scatter, plot_permutation,
        plot_volcano, plot_diff_bar,
        run_anova_kruskal, plot_multi_group_heatmap, plot_multi_group_boxplot,
        GROUP_COLORS, VOLCANO_COLORS,
    )
except ImportError:
    detect_species_columns_from_diff = build_abundance_matrix = compute_log2fc = None
    run_pca = run_plsda = run_oplsda = permutation_test_plsda = None
    plot_scores_scatter = plot_loading_scatter = plot_permutation = None
    plot_volcano = plot_diff_bar = None
    run_anova_kruskal = plot_multi_group_heatmap = plot_multi_group_boxplot = None
    GROUP_COLORS = VOLCANO_COLORS = None
import streamlit as st
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go

# ====================== KEGG 通路名映射（常用通路，不需要联网）======================
_PATHWAY_NAME_MAP = {
    'map01100': 'Metabolic pathways',
    'map01110': 'Biosynthesis of secondary metabolites',
    'map01120': 'Biosynthesis of secondary metabolites - bacteria',
    'map00260': 'Glycine, serine and threonine metabolism',
    'map00270': 'Cysteine and methionine metabolism',
    'map00280': 'Valine, leucine and isoleucine degradation',
    'map00290': 'Valine, leucine and isoleucine biosynthesis',
    'map00300': 'Lysine biosynthesis',
    'map00310': 'Lysine degradation',
    'map00330': 'Arginine biosynthesis',
    'map00340': 'Histidine metabolism',
    'map00350': 'Alanine, aspartate and glutamate metabolism',
    'map00360': 'Phenylalanine metabolism',
    'map00380': 'Tryptophan metabolism',
    'map00400': 'Phenylalanine, tyrosine and tryptophan biosynthesis',
    'map00430': 'Phenylpropanoid biosynthesis',
    'map00440': 'Phenylpropanoid biosynthesis - flavonoids',
    'map00450': 'Flavonoid biosynthesis',
    'map00460': 'Porphyrin metabolism',
    'map00480': 'Glutathione metabolism',
    'map00500': 'Starch and sucrose metabolism',
    'map00510': 'N-Glycan biosynthesis',
    'map00511': 'Other glycan degradation',
    'map00520': 'Amino sugar and nucleotide sugar metabolism',
    'map00522': 'Biosynthesis of unsaturated fatty acids',
    'map00523': 'alpha-Linolenic acid metabolism',
    'map00524': 'Phenylpropanoid biosynthesis - general',
    'map00530': 'Aminosugars metabolism',
    'map00531': 'Glycosaminoglycan degradation',
    'map00532': 'Glycosaminoglycan biosynthesis - chondroitin sulfate',
    'map00533': 'Glycosaminoglycan biosynthesis - keratan sulfate',
    'map00534': 'Glycosaminoglycan biosynthesis - heparan sulfate',
    'map00540': 'Lipopolysaccharide biosynthesis',
    'map00550': 'Peptidoglycan biosynthesis',
    'map00560': 'Pentose and glucuronate interconversions',
    'map00561': 'Glycerolipid metabolism',
    'map00562': 'Inositol phosphate metabolism',
    'map00563': 'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis',
    'map00564': 'Glycerophospholipid metabolism',
    'map00565': 'Ether lipid metabolism',
    'map00590': 'Arachidonic acid metabolism',
    'map00591': 'Linoleic acid metabolism',
    'map00592': 'alpha-Linolenic acid metabolism',
    'map00600': 'Sphingolipid metabolism',
    'map00620': 'Pentose phosphate pathway',
    'map00630': 'Glyoxylate and dicarboxylate metabolism',
    'map00640': 'Propanoate metabolism',
    'map00650': 'Butanoate metabolism',
    'map00660': 'C5-Branched dibasic acid metabolism',
    'map00670': 'One carbon pool by folate',
    'map00680': 'Methane metabolism',
    'map00700': 'Carbon fixation in photosynthetic organisms',
    'map00710': 'Carbon fixation pathways',
    'map00720': 'Carbon fixation in prokaryotes',
    'map00730': 'Thiamine metabolism',
    'map00740': 'Riboflavin metabolism',
    'map00750': 'Vitamin B6 metabolism',
    'map00760': 'Nicotinate and nicotinamide metabolism',
    'map00770': 'Pantothenate and CoA biosynthesis',
    'map00780': 'Biotin metabolism',
    'map00785': 'Lipoic acid metabolism',
    'map00790': 'Folate biosynthesis',
    'map00830': 'Selenocompound metabolism',
    'map00860': 'Porphyrin and chlorophyll metabolism',
    'map00900': 'Terpenoid backbone biosynthesis',
    'map00901': 'Camphor and pinene degradation',
    'map00902': 'Monoterpenoid biosynthesis',
    'map00903': 'Limonene and pinene degradation',
    'map00904': 'Diterpenoid biosynthesis',
    'map00905': 'Brassinosteroid biosynthesis',
    'map00906': 'Carotenoid biosynthesis',
    'map00908': 'Zeatin biosynthesis',
    'map00909': 'Sesquiterpenoid and triterpenoid biosynthesis',
    'map00910': 'Nitrogen metabolism',
    'map00920': 'Sulfur metabolism',
    'map00941': 'Flavonoid biosynthesis',
    'map00942': 'Anthocyanin biosynthesis',
    'map00943': 'Isoflavonoid biosynthesis',
    'map00944': 'Flavone and flavonol biosynthesis',
    'map00945': 'Anthocyanin biosynthesis (phenylpropanoid)',
    'map00950': 'Isoquinoline alkaloid biosynthesis',
    'map00960': 'Tropane, piperidine and pyridine alkaloid biosynthesis',
    'map00965': 'Betalain biosynthesis',
    'map00966': 'Glucosinolate biosynthesis',
    'map00970': 'Aminoacyl-tRNA biosynthesis',
    'map00980': 'Metabolism of xenobiotics by cytochrome P450',
    'map00981': 'Insect hormone biosynthesis',
    'map00982': 'Drug metabolism - cytochrome P450',
    'map00983': 'Drug metabolism - other enzymes',
    'map00984': 'Steroid degradation',
    'map00997': 'Other biosynthesis - alkaloids',
    'map00999': 'Biosynthesis of plant secondary metabolites',
    'map01040': 'Biosynthesis of unsaturated fatty acids',
    'map01060': 'Biosynthesis of plant secondary metabolites (aromatic)',
    'map01070': 'Biosynthesis of plant hormones',
    'map02010': 'ABC transporters',
    'map02020': 'Two-component system',
    'map02024': 'Quorum sensing',
    'map03010': 'Ribosome',
    'map03013': 'RNA transport',
    'map03015': 'mRNA surveillance pathway',
    'map03018': 'RNA degradation',
    'map03020': 'Ribosome biogenesis',
    'map03022': 'Basal transcription factors',
    'map03030': 'DNA replication',
    'map03040': 'Spliceosome',
    'map03060': 'Protein export',
    'map03008': 'Ribosome biogenesis in eukaryotes',
    'map03010': 'Ribosome',
    'map04070': 'Phosphatidylinositol signaling system',
    'map04075': 'Plant hormone signal transduction',
    'map04076': 'Plant-pathogen interaction',
    'map04120': 'Ubiquitin mediated proteolysis',
    'map04122': 'Sulfur relay system',
    'map04130': 'SNARE interactions in vesicular transport',
    'map04136': 'Autophagy - yeast',
    'map04140': 'Regulation of autophagy',
    'map04145': 'Phagosome',
    'map04150': 'mTOR signaling pathway',
    'map04210': 'Apoptosis',
    'map04260': 'Cardiac muscle contraction',
    'map04261': 'Adrenergic signaling in cardiomyocytes',
    'map04270': 'Vascular smooth muscle contraction',
    'map04310': 'Wnt signaling pathway',
    'map04330': 'Notch signaling pathway',
    'map04340': 'Hedgehog signaling pathway',
    'map04350': 'TGF-beta signaling pathway',
    'map04390': 'Hippo signaling pathway',
    'map04391': 'Hippo signaling pathway - fly',
    'map04510': 'Focal adhesion',
    'map04512': 'ECM-receptor interaction',
    'map04514': 'Cell adhesion molecules (CAMs)',
    'map04520': 'Adherens junction',
    'map04530': 'Tight junction',
    'map04540': 'Gap junction',
    'map04610': 'Complement and coagulation cascades',
    'map04611': 'Platelet activation',
    'map04612': 'Antigen processing and presentation',
    'map04614': 'Renin-angiotensin system',
    'map04620': 'Toll-like receptor signaling pathway',
    'map04621': 'NOD-like receptor signaling pathway',
    'map04622': 'RIG-I-like receptor signaling pathway',
    'map04623': 'Cytosolic DNA-sensing pathway',
    'map04625': 'C-type lectin receptor signaling pathway',
    'map04626': 'Plant-type R gene-mediated defense',
    'map04630': 'JAK-STAT signaling pathway',
    'map04640': 'Hematopoietic cell lineage',
    'map04650': 'Natural killer cell mediated cytotoxicity',
    'map04657': 'IL-17 signaling pathway',
    'map04658': 'Th1 and Th2 cell differentiation',
    'map04659': 'Th17 cell differentiation',
    'map04660': 'T cell receptor signaling pathway',
    'map04661': 'T cell receptor signaling pathway',
    'map04662': 'B cell receptor signaling pathway',
    'map04664': 'Fc epsilon RI signaling pathway',
    'map04666': 'Fc gamma R-mediated phagocytosis',
    'map04670': 'Leukocyte transendothelial migration',
    'map04672': 'Intestinal immune network for IgA production',
    'map04710': 'Circadian rhythm - fly',
    'map04713': 'Circadian entrainment',
    'map04714': 'Thermogenesis',
    'map04720': 'Long-term potentiation',
    'map04722': 'Neurotrophin signaling pathway',
    'map04723': 'Retrograde endocannabinoid signaling',
    'map04724': 'Glutamatergic synapse',
    'map04725': 'Cholinergic synapse',
    'map04726': 'Serotonergic synapse',
    'map04727': 'GABAergic synapse',
    'map04728': 'Dopaminergic synapse',
    'map04729': 'Long-term depression',
    'map04730': 'Long-term depression (cerebellum)',
    'map04740': 'Olfactory transduction',
    'map04742': 'Taste transduction',
    'map04744': 'Phototransduction',
    'map04750': 'Inflammatory mediator regulation of TRP channels',
    'map04810': 'Regulation of actin cytoskeleton',
    'map04814': 'Cytoskeleton in muscle cells',
    'map04910': 'Adipocytokine signaling pathway',
    'map04911': 'Insulin secretion',
    'map04912': 'GnRH signaling pathway',
    'map04913': 'Ovarian steroidogenesis',
    'map04914': 'Progesterone-mediated oocyte maturation',
    'map04915': 'Estrogen signaling pathway',
    'map04916': 'Melanogenesis',
    'map04917': 'Prolactin signaling pathway',
    'map04918': 'Thyroid hormone synthesis',
    'map04919': 'Thyroid hormone signaling pathway',
    'map04920': 'Adipocytokine signaling pathway',
    'map04921': 'Oxytocin signaling pathway',
    'map04922': 'Vasopressin-regulated water reabsorption',
    'map04923': 'Regulation of lipolysis in adipocytes',
    'map04924': 'Renin secretion',
    'map04925': 'Aldosterone synthesis and secretion',
    'map04926': 'Relaxin signaling pathway',
    'map04927': 'Cortisol synthesis and secretion',
    'map04928': 'Parathyroid hormone synthesis, secretion and action',
    'map04929': 'GnRH secretion',
    'map04930': 'Type II diabetes mellitus',
    'map04931': 'Insulin resistance',
    'map04932': 'Non-alcoholic fatty liver disease (NAFLD)',
    'map04933': 'AGE-RAGE signaling pathway in diabetic complications',
    'map04940': 'Type I diabetes mellitus',
    'map04950': 'Maturity onset diabetes of the young',
    'map04960': 'Aldosterone-regulated sodium reabsorption',
    'map04961': 'Endocrine and other factor-regulated calcium reabsorption',
    'map04962': 'Vasopressin-regulated water reabsorption',
    'map04964': 'Proximal tubule bicarbonate reclamation',
    'map04966': 'Collecting duct acid secretion',
    'map04970': 'Salivary secretion',
    'map04971': 'Gastric acid secretion',
    'map04972': 'Pancreatic secretion',
    'map04973': 'Carbohydrate digestion and absorption',
    'map04974': 'Protein digestion and absorption',
    'map04975': 'Fat digestion and absorption',
    'map04976': 'Bile secretion',
    'map04977': 'Vitamin digestion and absorption',
    'map04978': 'Mineral absorption',
    'map04979': 'Cholesterol metabolism',
    'map04980': 'Hemostasis',
    'map04981': 'Lipid digestion, absorption, transport, and metabolism',
    'map04982': 'Complement and coagulation cascades',
    'map04983': 'Drug metabolism - cytochrome P450',
    'map04984': 'Retinol metabolism',
    'map04985': 'Lactose synthesis and regulation',
    'map04986': 'Thermogenesis',
    'map04987': 'Fat digestion and absorption',
    'map04988': 'Linoleic acid metabolism',
    'map04989': 'alpha-Linolenic acid metabolism',
    'map04990': 'Phosphatidylinositol signaling system',
    'map04991': 'Inositol phosphate metabolism',
    'map04992': 'Sphingolipid metabolism',
    'map04993': 'Glycerolipid metabolism',
    'map04994': 'Glycerophospholipid metabolism',
    'map04995': 'Arachidonic acid metabolism',
    'map04996': 'Ether lipid metabolism',
    'map04997': 'Biosynthesis of unsaturated fatty acids',
    'map05010': 'Alzheimer disease',
    'map05012': 'Parkinson disease',
    'map05014': 'Amyotrophic lateral sclerosis (ALS)',
    'map05016': 'Huntington disease',
    'map05020': 'Prion diseases',
    'map05022': 'Pathways of neurodegeneration - multiple diseases',
    'map05030': 'Cocaine addiction',
    'map05031': 'Amphetamine addiction',
    'map05032': 'Morphine addiction',
    'map05033': 'Nicotine addiction',
    'map05034': 'Alcoholism',
    'map05010': 'Alzheimer disease',
    'map05100': 'Bacterial invasion of epithelial cells',
    'map05110': 'Vibrio cholerae infection',
    'map05111': 'Biofilm formation - Vibrio cholerae',
    'map05112': 'Pathogenic Escherichia coli infection',
    'map05113': 'Escherichia coli infection',
    'map05120': 'Epithelial cell signaling in Helicobacter pylori infection',
    'map05130': 'Pathogenic Escherichia coli infection',
    'map05131': 'Shigellosis',
    'map05132': 'Salmonella infection',
    'map05133': 'Pertussis',
    'map05134': 'Legionellosis',
    'map05135': 'Yersinia infection',
    'map05140': 'Leishmaniasis',
    'map05142': 'Chagas disease (American trypanosomiasis)',
    'map05143': 'African trypanosomiasis',
    'map05144': 'Malaria',
    'map05145': 'Toxoplasmosis',
    'map05146': 'Amoebiasis',
    'map05150': 'Staphylococcus aureus infection',
    'map05152': 'Tuberculosis',
    'map05160': 'Hepatitis C',
    'map05161': 'Hepatitis B',
    'map05162': 'Measles',
    'map05163': 'Human cytomegalovirus infection',
    'map05164': 'Influenza A',
    'map05165': 'Human papillomavirus infection',
    'map05166': 'Human T-cell leukemia virus type 1 infection',
    'map05167': 'Human papillomavirus type 16 infection',
    'map05168': 'Herpes simplex virus 1 infection',
    'map05169': 'Epstein-Barr virus infection',
    'map05170': 'Human immunodeficiency virus 1 infection',
    'map05171': 'Coronavirus disease - COVID-19',
    'map05175': 'Poliovirus infection',
    'map05176': 'Japanese encephalitis',
    'map05177': 'Hepatitis A',
    'map05178': 'Hepatitis D',
    'map05179': 'Hepatitis E',
    'map05200': 'Pathways in cancer',
    'map05202': 'Transcriptional misregulation in cancer',
    'map05203': 'Viral carcinogenesis',
    'map05204': 'Chemical carcinogenesis - activation',
    'map05205': 'Proteoglycans in cancer',
    'map05206': 'MicroRNAs in cancer',
    'map05207': 'Chemical carcinogenesis - receptor activation',
    'map05208': 'Chemical carcinogenesis - DNA adduct formation',
    'map05210': 'Colorectal cancer',
    'map05211': 'Renal cell carcinoma',
    'map05212': 'Pancreatic cancer',
    'map05213': 'Endometrial cancer',
    'map05214': 'Glioma',
    'map05215': 'Prostate cancer',
    'map05216': 'Thyroid cancer',
    'map05217': 'Basal cell carcinoma',
    'map05218': 'Melanoma',
    'map05219': 'Bladder cancer',
    'map05220': 'Chronic myeloid leukemia',
    'map05221': 'Acute myeloid leukemia',
    'map05222': 'Small cell lung cancer',
    'map05223': 'Non-small cell lung cancer',
    'map05224': 'Breast cancer',
    'map05225': 'Hepatocellular carcinoma',
    'map05226': 'Gastric cancer',
    'map05230': 'Central carbon metabolism in cancer',
    'map05231': 'Choline metabolism in cancer',
    'map05235': 'PD-L1 expression and PD-1 checkpoint pathway in cancer',
    'map05330': 'Allograft rejection',
    'map05331': 'Graft-versus-host disease',
    'map05332': 'Autoimmune thyroid disease',
    'map05340': 'Primary immunodeficiency',
    'map05350': 'Asthma',
    'map05360': 'Systemic lupus erythematosus',
    'map05410': 'Hypertrophic cardiomyopathy (HCM)',
    'map05412': 'Arrhythmogenic right ventricular cardiomyopathy (ARVC)',
    'map05414': 'Dilated cardiomyopathy (DCM)',
    'map05416': 'Viral myocarditis',
    'map05418': 'Fluid shear stress and atherosclerosis',
    'map01200': 'Carbon metabolism',
    'map01210': '2-Oxocarboxylic acid metabolism',
    'map01212': 'Fatty acid metabolism',
    'map01230': 'Biosynthesis of amino acids',
    'map01232': 'Nucleotide metabolism',
    'map01240': 'Biosynthesis of cofactors',
    'map01250': 'Biosynthesis of alkaloids',
    'map01220': 'Degradation of aromatic compounds',
}
st.set_page_config(
    page_title="差异代谢物药用筛选",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ====== MetaboLab Aura Design System CSS ======
# 参考：stitch_high_end_website_prototype - Tailwind MD3 原型设计精确复刻

st.markdown("""
<style>
/* ====== Google Fonts: Inter ====== */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800;900&family=Noto+Sans+SC:wght@400;500;700&display=swap');

/* ====== 全局基础 ====== */
html, body { font-family: "Inter", "Noto Sans SC", -apple-system, sans-serif !important; -webkit-font-smoothing: antialiased; letter-spacing: -0.01em; }
h1, h2, h3, h4 { letter-spacing: -0.02em !important; }

/* ====== 根背景 - 薰衣草白 ====== */
.stApp { background-color: #fef7ff !important; }

/* ====== 顶部导航栏 - 玻璃拟态 ====== */
[data-testid="stHeader"] { background: rgba(254, 247, 255, 0.82) !important; backdrop-filter: blur(24px) !important; -webkit-backdrop-filter: blur(24px) !important; border-bottom: none !important; box-shadow: 0 40px 60px -15px rgba(99, 14, 212, 0.05) !important; }

/* ====== 侧边栏 - 薰衣草紫 ====== */
[data-testid="stSidebar"] { background-color: #f3ebfa !important; border-right: none !important; }
[data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3 { color: #630ed4 !important; font-weight: 700 !important; }
[data-testid="stSidebar"] hr { border: none !important; border-top: 1px solid rgba(204,195,216,0.4) !important; }
[data-testid="stSidebar"] .stMarkdown p, [data-testid="stSidebar"] span, [data-testid="stSidebar"] label { color: #1d1a24 !important; }

/* ====== 侧边栏 Selectbox ====== */
[data-testid="stSidebar"] [data-baseweb="select"] > div, [data-testid="stSidebar"] [data-baseweb="multiselect"] > div { background-color: #ffffff !important; border: none !important; border-radius: 1rem !important; box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important; }
[data-testid="stSidebar"] .stSelectbox > label, [data-testid="stSidebar"] .stMultiSelect > label, [data-testid="stSidebar"] label { color: #4a4455 !important; font-weight: 600 !important; font-size: 11px !important; letter-spacing: 0.06em !important; text-transform: uppercase !important; }

/* ====== 侧边栏 Slider ====== */
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-track { background: linear-gradient(90deg, #630ed4, #7c3aed) !important; height: 4px !important; border-radius: 9999px !important; border: none !important; }
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-thumb { background-color: #630ed4 !important; width: 16px !important; height: 16px !important; box-shadow: 0 2px 8px rgba(99,14,212,0.4) !important; }
[data-testid="stSidebar"] [data-baseweb="slider"] .MuiSlider-rail { background-color: #e8dfee !important; height: 4px !important; border-radius: 9999px !important; opacity: 1 !important; }

/* ====== Tab 胶囊式 ====== */
.stTabs [data-baseweb="tab-list"] { background-color: #f3ebfa !important; border-radius: 1.5rem !important; padding: 6px !important; gap: 4px !important; border: none !important; }
.stTabs [data-baseweb="tab"] { border-radius: 1rem !important; font-weight: 500 !important; font-size: 13px !important; color: #4a4455 !important; background: transparent !important; border: none !important; padding: 8px 18px !important; transition: all 0.25s ease !important; }
.stTabs [data-baseweb="tab"]:hover { background-color: rgba(99,14,212,0.08) !important; color: #630ed4 !important; }
.stTabs [aria-selected="true"] { background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important; color: #ffffff !important; font-weight: 600 !important; box-shadow: 0 4px 20px rgba(99,14,212,0.4) !important; border: none !important; }

/* ====== 主按钮 ====== */
.stButton > button[kind="primary"], div[data-testid="stMainBlockContainer"] button[kind="primary"] { background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important; color: #ffffff !important; border: none !important; border-radius: 9999px !important; padding: 0.6rem 2rem !important; font-weight: 600 !important; font-size: 14px !important; box-shadow: 0 8px 16px -4px rgba(99,14,212,0.35) !important; transition: all 0.2s ease !important; }
.stButton > button[kind="primary"]:hover { background: linear-gradient(135deg, #5509c6 0%, #6d28d9 100%) !important; box-shadow: 0 12px 24px -4px rgba(99,14,212,0.45) !important; transform: translateY(-1px) !important; opacity: 0.95 !important; }
.stButton > button[kind="primary"]:active { transform: scale(0.98) translateY(0) !important; }

/* ====== 普通按钮 ====== */
.stButton > button { border-radius: 1rem !important; font-weight: 500 !important; font-size: 13px !important; border: none !important; background-color: #f3ebfa !important; color: #1d1a24 !important; padding: 0.5rem 1.2rem !important; transition: all 0.2s ease !important; }
.stButton > button:hover { background-color: #ede5f4 !important; box-shadow: 0 2px 12px rgba(99,14,212,0.12) !important; }

/* ====== 下载按钮 ====== */
.stDownloadButton > button { background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important; color: #ffffff !important; border: none !important; border-radius: 9999px !important; font-weight: 600 !important; font-size: 13px !important; padding: 0.5rem 1.5rem !important; box-shadow: 0 4px 16px rgba(99,14,212,0.3) !important; }
.stDownloadButton > button:hover { box-shadow: 0 6px 24px rgba(99,14,212,0.4) !important; transform: translateY(-1px) !important; }

/* ====== Metric ====== */
[data-testid="stMetricValue"] { color: #630ed4 !important; font-weight: 800 !important; font-size: 36px !important; letter-spacing: -0.04em !important; }
[data-testid="stMetricLabel"] { color: #4a4455 !important; font-weight: 500 !important; font-size: 12px !important; letter-spacing: 0.03em !important; text-transform: uppercase !important; }

/* ====== DataFrame ====== */
.stDataFrame { border-radius: 1.5rem !important; overflow: hidden !important; border: none !important; box-shadow: 0 20px 40px -15px rgba(99,14,212,0.06) !important; }
.stDataFrame thead tr th { background-color: #f3ebfa !important; color: #4a4455 !important; font-weight: 600 !important; font-size: 11px !important; letter-spacing: 0.08em !important; text-transform: uppercase !important; border-bottom: none !important; padding: 14px 20px !important; }
.stDataFrame thead th:first-child { border-radius: 1.5rem 0 0 0 !important; }
.stDataFrame thead th:last-child { border-radius: 0 1.5rem 0 0 !important; }
.stDataFrame tbody tr { background-color: #ffffff !important; transition: background-color 0.15s ease !important; }
.stDataFrame tbody tr:hover { background-color: rgba(249,241,255,0.6) !important; }
.stDataFrame tbody tr:nth-child(even) { background-color: rgba(237,229,244,0.15) !important; }
.stDataFrame tbody tr:nth-child(even):hover { background-color: rgba(249,241,255,0.6) !important; }
.stDataFrame tbody td { color: #1d1a24 !important; font-size: 13px !important; padding: 12px 20px !important; border-bottom: none !important; }

/* ====== 展开面板 ====== */
details { background-color: #ffffff !important; border-radius: 1.5rem !important; border: none !important; box-shadow: 0 4px 20px rgba(99,14,212,0.07) !important; overflow: hidden !important; }
details summary { background-color: #f9f1ff !important; color: #1d1a24 !important; font-weight: 600 !important; font-size: 13px !important; border-radius: 1.5rem !important; padding: 14px 20px !important; }
details[open] summary { border-radius: 1.5rem 1.5rem 0 0 !important; }
details > div { background-color: #ffffff !important; padding: 16px 20px !important; }

/* ====== Alert ====== */
.stAlert { border-radius: 1.5rem !important; border: none !important; box-shadow: 0 4px 20px rgba(99,14,212,0.08) !important; }

/* ====== Selectbox / Multiselect ====== */
.stSelectbox [data-baseweb="select"] > div, .stMultiSelect [data-baseweb="select"] > div { background-color: #ffffff !important; border: none !important; border-radius: 1rem !important; box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important; }
.stSelectbox label, .stMultiSelect label { color: #4a4455 !important; font-weight: 600 !important; font-size: 11px !important; letter-spacing: 0.06em !important; text-transform: uppercase !important; }

/* ====== Slider ====== */
.stSlider [data-baseweb="slider"] .MuiSlider-track { background: linear-gradient(90deg, #630ed4, #7c3aed) !important; height: 4px !important; border-radius: 9999px !important; border: none !important; }
.stSlider [data-baseweb="slider"] .MuiSlider-thumb { background-color: #630ed4 !important; width: 18px !important; height: 18px !important; box-shadow: 0 2px 10px rgba(99,14,212,0.45) !important; }
.stSlider [data-baseweb="slider"] .MuiSlider-rail { background-color: #e8dfee !important; height: 4px !important; border-radius: 9999px !important; opacity: 1 !important; }

/* ====== Text / Number Input ====== */
.stTextInput input, .stTextArea textarea, .stNumberInput input { background-color: #ffffff !important; border: none !important; border-radius: 1rem !important; box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important; font-size: 13px !important; color: #1d1a24 !important; padding: 8px 14px !important; }
.stTextInput input:focus, .stTextArea textarea:focus, .stNumberInput input:focus { outline: 2px solid #7c3aed !important; box-shadow: 0 2px 12px rgba(99,14,212,0.15) !important; }

/* ====== Spinner ====== */
.stSpinner > div { border-color: rgba(99,14,212,0.2) !important; border-top-color: #630ed4 !important; }

/* ====== Progress ====== */
.stProgress > div > div > div { background: linear-gradient(90deg, #630ed4, #7c3aed) !important; border-radius: 9999px !important; height: 6px !important; }

/* ====== Plotly ====== */
.js-plotly-plot .plotly { border-radius: 1.5rem !important; overflow: hidden !important; box-shadow: 0 4px 24px rgba(99,14,212,0.1) !important; }
.js-plotly-plot .plotly .modebar { background: rgba(254,247,255,0.9) !important; border-radius: 1rem !important; padding: 4px 8px !important; }

/* ====== 全局滚动条 ====== */
::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb { background: #ccc3d8; border-radius: 9999px; }
::-webkit-scrollbar-thumb:hover { background: #7b7487; }

/* ====== Container 卡片 ====== */
[data-testid="stVerticalBlock"] [data-testid="stHorizontalBlock"] > div { background-color: #ffffff !important; border-radius: 1.5rem !important; box-shadow: 0 2px 16px rgba(99,14,212,0.06) !important; padding: 20px !important; }

/* ====== 文件上传器 ====== */
[data-testid="stFileUploader"] > div { border-radius: 1.5rem !important; border: 2px dashed rgba(99,14,212,0.25) !important; background-color: rgba(249,241,255,0.5) !important; backdrop-filter: blur(8px) !important; transition: all 0.3s ease !important; }
[data-testid="stFileUploader"] > div:hover { border-color: rgba(99,14,212,0.5) !important; background-color: rgba(249,241,255,0.8) !important; }
</style>
""", unsafe_allow_html=True)

plt.rcParams['font.sans-serif'] = ['SimHei', 'WenQuanYi Micro Hei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ====================== 药理分类规则 ======================
PHARMA_SUPER_CLASS_KEYWORDS = [
    'Alkaloids', 'Phenylpropanoids', 'Lignans', 'Benzenoids',
    'Organoheterocyclic', 'Organosulfur', 'Organonitrogen',
    'Organic acids', 'Nucleosides',
]

PHARMA_COMPOUND_KEYWORDS = [
    'Peptides', 'Steroids', 'Vitamins', 'Antibiotics',
]

PHARMA_KEGG_KEYWORDS = [
    'Amino acid metabolism',
    'Biosynthesis of other secondary metabolites',
    'Biosynthesis of alkaloids',
    'Biosynthesis of phenylpropanoids',
    'Biosynthesis of terpenoids and polyketides',
    'Metabolism of terpenoids and polyketides',
    'Biosynthesis of antibiotics',
    'Biosynthesis of amino acids',
]

# 明星分子药理评分
PHARMA_SCORES = {
    'Alkaloids and derivatives': 1.0,
    'Alkaloids': 1.0,
    'Flavonoids': 0.9,
    'Terpenoids': 0.8,
    'Lignans, neolignans and related compounds': 0.8,
    'Phenylpropanoids and polyketides': 0.7,
    'Benzenoids': 0.6,
    'Organoheterocyclic compounds': 0.5,
    'Organic acids and derivatives': 0.4,
    'Nucleosides, nucleotides, and analogues': 0.4,
    'Organic oxygen compounds': 0.3,
    'Lipids and lipid-like molecules': 0.2,
}

# ====================== 核心处理函数 ======================

def read_uploaded_file(uploaded_file) -> pd.DataFrame:
    """根据文件扩展名自动推断格式，读取上传的文件为 DataFrame"""
    name = uploaded_file.name
    content = uploaded_file.read()

    # 判断是否为 TSV（检查第一行是否含制表符）
    first_line = content.split(b'\n')[0]
    if b'\t' in first_line:
        df = pd.read_csv(io.BytesIO(content), sep='\t')
    elif name.endswith('.xlsx') or name.endswith('.xls'):
        try:
            df = pd.read_excel(io.BytesIO(content), engine='xlrd')
        except Exception:
            df = pd.read_excel(io.BytesIO(content), engine='openpyxl')
    else:
        df = pd.read_csv(io.BytesIO(content), sep='\t')

    df.columns = df.columns.str.strip()
    return df


def detect_species_columns(df) -> dict:
    """扫描 DataFrame 列名，检测各物种的丰度列"""
    all_cols = df.columns.tolist()
    species_keywords = {
        'LJ': 'LJ_Root', 'BC': 'BC_Root', 'XY': 'XY_Root',
        'EH': 'EH_Root', 'MB': 'MB_Root',
    }
    result = {sp: [] for sp in species_keywords}
    for col in all_cols:
        for sp, kw in species_keywords.items():
            if kw in col:
                result[sp].append(col)
    return result


def process_diff_files(dfs: list) -> pd.DataFrame:
    """合并多个 diff DataFrame，按 VIP 降序去重"""
    for df in dfs:
        df.columns = df.columns.str.strip()
    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.sort_values('Vip_plsda', ascending=False)
    merged = merged.drop_duplicates(subset=['Metabolite'], keep='first')
    merged = merged.reset_index(drop=True)
    return merged


def filter_differential(df, vip_thresh=1.0, p_thresh=0.05):
    """VIP > threshold AND P_value < threshold"""
    mask = (df['Vip_plsda'] > vip_thresh) & (df['P_value'] < p_thresh)
    return df.loc[mask].copy().reset_index(drop=True)


def annotate_metabolites(diff_df, anno_df):
    """将 anno 注释关联到 diff 数据，按 Metabolite ~ metab 匹配"""
    annotated = diff_df.merge(
        anno_df,
        left_on='Metabolite',
        right_on='metab',
        how='left'
    )
    return annotated


def is_pharmacological(row) -> tuple:
    """判断化合物是否具有潜在药理活性，返回 (bool, reason)"""
    sc = str(row.get('super_class', '-')).strip()
    cc = str(row.get('compound_first_category', '-')).strip()
    kegg2 = str(row.get('kegg_second_category', '-')).strip()
    kegg1 = str(row.get('kegg_first_category', '-')).strip()
    name = str(row.get('compound_name', '-')).strip()

    if name not in ('-', '', 'nan'):
        for kw in PHARMA_SUPER_CLASS_KEYWORDS:
            if kw.lower() in sc.lower():
                return True, f"super_class={kw}"
        for kw in PHARMA_COMPOUND_KEYWORDS:
            if kw.lower() in cc.lower():
                return True, f"compound_category={kw}"
        for kw in PHARMA_KEGG_KEYWORDS:
            if kw.lower() in kegg2.lower():
                return True, f"kegg_pathway={kw}"
        if 'Biosynthesis of secondary metabolites' in kegg1:
            return True, "Biosynthesis of secondary metabolites"
        if 'Metabolism of terpenoids and polyketides' in kegg1:
            return True, "terpenoids/polyketides metabolism"

    if name in ('-', '', 'nan'):
        if sc not in ('-', '', 'nan', 'Not Available'):
            if any(kw.lower() in sc.lower() for kw in PHARMA_SUPER_CLASS_KEYWORDS):
                return True, f"super_class_only={sc}"

    return False, ''


def filter_pharmacological(df) -> pd.DataFrame:
    """筛选药理活性化合物"""
    results = []
    for _, row in df.iterrows():
        is_pharma, reason = is_pharmacological(row)
        if is_pharma:
            row_dict = row.to_dict()
            row_dict['_pharma_reason'] = reason
            results.append(row_dict)
    return pd.DataFrame(results)


def calculate_abundance(pharma_metabolites: list, diff_df: pd.DataFrame) -> pd.DataFrame:
    """
    计算每个候选化合物的各物种平均丰度。
    diff_df: 合并去重后的完整 diff DataFrame
    """
    species_cols = detect_species_columns(diff_df)
    all_cols = diff_df.columns.tolist()

    rows = []
    for metab in pharma_metabolites:
        row_data = {'Metabolite': metab}
        metab_row = diff_df[diff_df['Metabolite'] == metab]

        if metab_row.empty:
            for sp in species_cols:
                row_data[f'{sp}_mean_abundance'] = None
            rows.append(row_data)
            continue

        for sp, cols in species_cols.items():
            if cols:
                vals = []
                for col in cols:
                    if col in all_cols:
                        v = metab_row[col].values[0]
                        if pd.notna(v):
                            try:
                                vals.append(float(v))
                            except (ValueError, TypeError):
                                pass
                if vals:
                    row_data[f'{sp}_mean_abundance'] = round(np.mean(vals), 6)
                else:
                    row_data[f'{sp}_mean_abundance'] = None
            else:
                row_data[f'{sp}_mean_abundance'] = None
        rows.append(row_data)

    return pd.DataFrame(rows)


def compute_star_scores(pharma_df: pd.DataFrame, abundance_df: pd.DataFrame) -> pd.DataFrame:
    """
    计算明星分子综合评分：
    综合得分 = VIP_score*0.3 + species_diff_score*0.3 + pharma_score*0.4
    """
    df = pharma_df.merge(abundance_df, on='Metabolite', how='left')

    vip_max = df['Vip_plsda'].max()
    df['_vip_norm'] = df['Vip_plsda'] / vip_max if vip_max > 0 else 0

    # 物种差异倍数（max / min，log处理）
    species_cols = [c for c in df.columns if c.endswith('_mean_abundance')]
    df['_max_abund'] = df[species_cols].replace(0, np.nan).max(axis=1)
    df['_min_abund'] = df[species_cols].replace(0, np.nan).min(axis=1)
    df['_diff_ratio'] = df['_max_abund'] / (df['_min_abund'] + 1e-9)
    df['_diff_ratio'] = df['_diff_ratio'].replace([np.inf, -np.inf], np.nan)
    diff_max = df['_diff_ratio'].quantile(0.95)
    df['_diff_norm'] = np.clip(df['_diff_ratio'] / diff_max, 0, 1)

    # 药理评分
    df['_pharma_score'] = df['super_class'].map(
        lambda x: PHARMA_SCORES.get(str(x).strip(), 0.3)
    )

    df['_star_score'] = (
        df['_vip_norm'] * 0.3 +
        df['_diff_norm'].fillna(0) * 0.3 +
        df['_pharma_score'] * 0.4
    )

    return df.sort_values('_star_score', ascending=False).reset_index(drop=True)


def extract_kegg_ids(df) -> list:
    """从 pathway_id 列提取所有 KEGG ID"""
    kegg_ids = []
    for _, row in df.iterrows():
        pid = str(row.get('pathway_id', '-'))
        for p in pid.split(';'):
            p = p.strip()
            if p.startswith('map') or p.startswith('ko'):
                kegg_ids.append(p)
    return list(set(kegg_ids))


# ====================== Excel 下载辅助 ======================

def df_to_excel_bytes(df: pd.DataFrame) -> bytes:
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Sheet1')
    output.seek(0)
    return output.getvalue()


# ====================== 可视化函数 ======================

def plot_heatmap(pharma_df: pd.DataFrame, abundance_df: pd.DataFrame, top_n=50):
    """绘制候选化合物丰度聚类热图"""
    from scipy import stats

    species_order = ['LJ', 'BC', 'XY', 'EH', 'MB']
    abundance_cols = [f'{sp}_mean_abundance' for sp in species_order
                      if f'{sp}_mean_abundance' in abundance_df.columns]

    if not abundance_cols:
        return None

    # 取 Top N VIP
    top_metabs = pharma_df.nlargest(top_n, 'Vip_plsda')['Metabolite'].tolist()
    plot_df = abundance_df[abundance_df['Metabolite'].isin(top_metabs)].copy()
    plot_df = plot_df.set_index('Metabolite').loc[top_metabs]

    # z-score 标准化
    z_df = plot_df[abundance_cols].apply(
        lambda x: stats.zscore(x, nan_policy='omit'), axis=0
    )

    # 简化列名
    z_df.columns = [c.replace('_mean_abundance', '') for c in z_df.columns]

    # 简化行名（太长则截断）
    z_df.index = [name[:30] + '...' if len(name) > 30 else name
                  for name in z_df.index]

    fig, ax = plt.subplots(figsize=(8, max(10, top_n * 0.25)))
    sns.heatmap(
        z_df, cmap='RdBu_r', center=0,
        xticklabels=True, yticklabels=True,
        cbar_kws={'label': 'Z-score'},
        ax=ax
    )
    ax.set_title(f'Candidate Compounds Heatmap (Top {top_n} by VIP)', fontsize=13)
    ax.set_xlabel('Species')
    ax.set_ylabel('Metabolite')
    plt.xticks(rotation=0)
    plt.yticks(fontsize=7)
    plt.tight_layout()
    return fig


def plot_boxplot(abundance_df: pd.DataFrame, pharma_df: pd.DataFrame, top_n=20):
    """绘制各物种丰度箱线图（宽转长）"""
    species_order = ['LJ', 'BC', 'XY', 'EH', 'MB']
    abundance_cols = [f'{sp}_mean_abundance' for sp in species_order
                      if f'{sp}_mean_abundance' in abundance_df.columns]

    if not abundance_cols:
        return None

    top_metabs = pharma_df.nlargest(top_n, 'Vip_plsda')['Metabolite'].tolist()
    plot_df = abundance_df[abundance_df['Metabolite'].isin(top_metabs)].copy()

    # 宽转长
    melt_df = plot_df.melt(
        id_vars=['Metabolite'],
        value_vars=abundance_cols,
        var_name='Species',
        value_name='Abundance'
    )
    melt_df['Species'] = melt_df['Species'].str.replace('_mean_abundance', '')
    melt_df['Abundance'] = pd.to_numeric(melt_df['Abundance'], errors='coerce')
    melt_df = melt_df.dropna(subset=['Abundance'])

    # 简化化合物名称
    melt_df['Metabolite_short'] = melt_df['Metabolite'].apply(
        lambda x: x[:25] + '...' if len(x) > 25 else x
    )

    fig = px.box(
        melt_df,
        x='Species',
        y='Abundance',
        color='Species',
        title=f'Abundance Distribution per Species (Top {top_n} VIP)',
        points='outliers',
        hover_data={'Metabolite': True, 'Abundance': ':.4f'}
    )
    fig.update_layout(
        boxmode='group',
        showlegend=True,
        xaxis_title='Species',
        yaxis_title='Mean Abundance (log2)',
        font=dict(size=12),
    )
    return fig


def plot_radar(pharma_row: dict, abundance_row: dict) -> go.Figure:
    """绘制单个化合物的五物种丰度雷达图"""
    species_order = ['LJ', 'BC', 'XY', 'EH', 'MB']
    abundance_cols = [f'{sp}_mean_abundance' for sp in species_order
                      if f'{sp}_mean_abundance' in abundance_row]

    vals = []
    labels = []
    for col in abundance_cols:
        sp = col.replace('_mean_abundance', '')
        v = abundance_row.get(col)
        if v is not None and not pd.isna(v):
            vals.append(float(v))
            labels.append(sp)
        else:
            vals.append(0)
            labels.append(sp)

    if len(vals) < 3:
        return None

    # 闭合雷达图
    vals = vals + [vals[0]]
    labels = labels + [labels[0]]

    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(
        r=vals,
        theta=labels,
        fill='toself',
        fillcolor='rgba(31,119,180,0.3)',
        line_color='rgb(31,119,180)',
        name=pharma_row.get('Metabolite', '')[:20],
    ))
    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True)),
        showlegend=False,
        title=dict(
            text=pharma_row.get('Metabolite', '')[:40],
            font=dict(size=12)
        ),
        height=300,
    )
    return fig

# ====================== Streamlit UI ======================

def main():
    # ===== 品牌头部横幅 =====
    st.markdown("""
    <div style="background: linear-gradient(135deg, #630ed4 0%, #7c3aed 50%, #9f67f5 100%);
                padding: 24px 32px; border-radius: 1.5rem; margin-bottom: 24px;
                box-shadow: 0 8px 32px rgba(99, 14, 212, 0.3);">
        <div style="display: flex; align-items: center; justify-content: space-between;">
            <div style="display: flex; align-items: center; gap: 16px;">
                <div style="font-size: 36px; filter: drop-shadow(0 2px 8px rgba(0,0,0,0.2));">&#129704;</div>
                <div>
                    <h1 style="color: #ffffff; margin: 0; font-size: 22px; font-weight: 700; letter-spacing: -0.03em; text-shadow: 0 1px 4px rgba(0,0,0,0.15);">
                        &#24046;&#24322;&#20195;&#35874;&#29289;&#33647;&#29992;&#31579;&#36873;&#24179;&#21488;
                    </h1>
                    <p style="color: rgba(255,255,255,0.82); margin: 6px 0 0 0; font-size: 13px; letter-spacing: 0.01em;">
                        Multivariate Stats &middot; Volcano Plot &middot; KEGG Enrichment &middot; Network Pharmacology
                    </p>
                </div>
            </div>
            <div style="background: rgba(255,255,255,0.15); border-radius: 1rem; padding: 8px 16px; backdrop-filter: blur(8px);">
                <span style="color: rgba(255,255,255,0.9); font-size: 12px; font-weight: 600; letter-spacing: 0.05em; text-transform: uppercase;">MetaboLab</span>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()

