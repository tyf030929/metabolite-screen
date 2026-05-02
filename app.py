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
    from network_pharma import (
        query_smiles_by_name, batch_query_smiles,
        query_swiss_target_prediction, batch_swiss_target_prediction,
        run_gseapy_enrichment, merge_enrichment_results, plot_enrichment_dotplot,
    )
except ImportError:
    query_smiles_by_name = batch_query_smiles = None
    query_swiss_target_prediction = batch_swiss_target_prediction = None
    run_gseapy_enrichment = merge_enrichment_results = plot_enrichment_dotplot = None

try:
    from pharma_cache import match_pharma_online, query_pharma_info, compute_pharma_evidence_score as compute_pharma_evidence_score_online, get_cache_status
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

# ====================== Aura Scientific 设计系统 CSS ======================
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;900&display=swap');
@import url('https://fonts.googleapis.com/css2?family=Noto+Sans+SC:wght@400;500;700&display=swap');

html, body, [class*="css-"] {
    font-family: 'Inter', 'Noto Sans SC', -apple-system, sans-serif !important;
    letter-spacing: -0.01em;
    -webkit-font-smoothing: antialiased;
}

h1, h2, h3, h4 { letter-spacing: -0.02em !important; }

.stApp { background-color: #fef7ff !important; }

[data-testid="stSidebar"] { background-color: #f3ebfa !important; border-right: none !important; }
[data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2,
[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] h1,
[data-testid="stSidebar"] .stTitle { color: #630ed4 !important; font-weight: 700 !important; }

.stTabs [data-baseweb="tab-list"] {
    background-color: #f3ebfa !important; border-radius: 1.5rem !important;
    padding: 6px !important; gap: 4px !important; border: none !important; box-shadow: none !important;
}
.stTabs [data-baseweb="tab"] {
    border-radius: 1rem !important; font-weight: 500 !important; font-size: 13px !important;
    color: #4a4455 !important; background: transparent !important; border: none !important;
    padding: 8px 16px !important; transition: all 0.2s ease !important;
}
.stTabs [data-baseweb="tab"]:hover { background-color: rgba(99,14,212,0.08) !important; color: #630ed4 !important; }
.stTabs [aria-selected="true"] {
    background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important;
    color: #ffffff !important; font-weight: 600 !important;
    box-shadow: 0 4px 20px rgba(99,14,212,0.35) !important;
}

.stButton > button[kind="primary"],
div[data-testid="stMainBlockContainer"] button[kind="primary"] {
    background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important;
    color: #ffffff !important; border: none !important; border-radius: 9999px !important;
    padding: 0.55rem 1.8rem !important; font-weight: 600 !important; font-size: 14px !important;
    box-shadow: 0 4px 20px rgba(99,14,212,0.35) !important; transition: all 0.2s ease !important;
}
.stButton > button[kind="primary"]:hover {
    background: linear-gradient(135deg, #5509c6 0%, #6d28d9 100%) !important;
    box-shadow: 0 6px 28px rgba(99,14,212,0.45) !important; transform: translateY(-1px) !important;
}

.stButton > button {
    border-radius: 1rem !important; font-weight: 500 !important; font-size: 13px !important;
    border: none !important; background-color: #f3ebfa !important; color: #1d1a24 !important;
    padding: 0.5rem 1.2rem !important; transition: all 0.2s ease !important; box-shadow: none !important;
}
.stButton > button:hover { background-color: #ede5f4 !important; }

.stDownloadButton > button {
    border-radius: 9999px !important; border: none !important;
    background: linear-gradient(135deg, #630ed4 0%, #7c3aed 100%) !important;
    color: #ffffff !important; font-weight: 600 !important; font-size: 13px !important;
    padding: 0.5rem 1.5rem !important; box-shadow: 0 4px 16px rgba(99,14,212,0.3) !important;
}
.stDownloadButton > button:hover { box-shadow: 0 6px 24px rgba(99,14,212,0.4) !important; transform: translateY(-1px) !important; }

[data-testid="stMetricValue"] { color: #630ed4 !important; font-weight: 700 !important; font-size: 28px !important; letter-spacing: -0.03em !important; }
[data-testid="stMetricLabel"] { color: #4a4455 !important; font-weight: 500 !important; font-size: 12px !important; }

.stDataFrame { border-radius: 1.5rem !important; overflow: hidden !important; border: none !important; box-shadow: 0 2px 16px rgba(99,14,212,0.08) !important; }
.stDataFrame thead th { background-color: #f3ebfa !important; color: #4a4455 !important; font-weight: 600 !important; font-size: 11px !important; letter-spacing: 0.08em !important; text-transform: uppercase !important; border-bottom: none !important; padding: 12px 16px !important; }
.stDataFrame tbody tr { background-color: #ffffff !important; transition: background-color 0.15s ease !important; }
.stDataFrame tbody tr:hover { background-color: #f9f1ff !important; }
.stDataFrame tbody td { color: #1d1a24 !important; font-size: 13px !important; padding: 10px 16px !important; border-bottom: 1px solid #f3ebfa !important; }
.stDataFrame thead th:first-child { border-radius: 1.5rem 0 0 0 !important; }
.stDataFrame thead th:last-child { border-radius: 0 1.5rem 0 0 !important; }

hr { border: none !important; margin: 0 !important; }

details { background-color: #ffffff !important; border-radius: 1.5rem !important; border: none !important; box-shadow: 0 2px 16px rgba(99,14,212,0.07) !important; overflow: hidden !important; }
details summary { background-color: #f3ebfa !important; color: #1d1a24 !important; font-weight: 600 !important; border-radius: 1.5rem !important; padding: 12px 20px !important; }
details[open] summary { border-radius: 1.5rem 1.5rem 0 0 !important; }
details > div { background-color: #ffffff !important; padding: 16px 20px !important; }

.stSelectbox [data-baseweb="select"] > div, .stMultiSelect [data-baseweb="select"] > div {
    background-color: #ffffff !important; border: none !important; border-radius: 1rem !important;
    box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important;
}
.stSelectbox label, .stMultiSelect label {
    color: #4a4455 !important; font-weight: 600 !important; font-size: 11px !important;
    letter-spacing: 0.05em !important; text-transform: uppercase !important;
}

.stSlider [data-baseweb="slider"] .MuiSlider-track {
    background: linear-gradient(90deg, #630ed4, #7c3aed) !important; border: none !important;
}
.stSlider [data-baseweb="slider"] .MuiSlider-thumb { background-color: #630ed4 !important; box-shadow: 0 2px 8px rgba(99,14,212,0.4) !important; }
.stSlider [data-baseweb="slider"] .MuiSlider-rail { background-color: #e8dfee !important; }

.stTextInput input, .stTextArea textarea {
    background-color: #ffffff !important; border: none !important; border-radius: 1rem !important;
    box-shadow: 0 1px 8px rgba(99,14,212,0.1) !important; font-size: 13px !important; color: #1d1a24 !important;
}
.stTextInput input:focus, .stTextArea textarea:focus {
    outline: 2px solid #7c3aed !important; box-shadow: 0 2px 12px rgba(99,14,212,0.2) !important;
}

.stProgress > div > div > div { background: linear-gradient(90deg, #630ed4, #7c3aed) !important; border-radius: 9999px !important; }

[data-testid="stVerticalBlock"] > div:has(> .stCard), .element-container > div > .stCard {
    background-color: #ffffff !important; border-radius: 1.5rem !important; border: none !important;
    box-shadow: 0 2px 16px rgba(99,14,212,0.07) !important; padding: 20px !important;
}

[data-testid="stImage"] figure { border-radius: 1.5rem !important; overflow: hidden !important; box-shadow: 0 4px 24px rgba(99,14,212,0.1) !important; }
.js-plotly-plot .plotly .modebar { background: rgba(254,247,255,0.9) !important; border-radius: 1rem !important; }

[data-testid="stVerticalBlock"] { scrollbar-width: thin; scrollbar-color: #e8dfee transparent; }
.stMain > div:first-child { background: linear-gradient(to bottom, rgba(99,14,212,0.02) 0%, transparent 300px) !important; }
[data-testid="stMarkdownContainer"] h1 { color: #ffffff !important; }
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
PHARMA_COMPOUND_KEYWORDS = ['Peptides', 'Steroids', 'Vitamins', 'Antibiotics']
PHARMA_KEGG_KEYWORDS = [
    'Amino acid metabolism', 'Biosynthesis of other secondary metabolites',
    'Biosynthesis of alkaloids', 'Biosynthesis of phenylpropanoids',
    'Biosynthesis of terpenoids and polyketides', 'Metabolism of terpenoids and polyketides',
    'Biosynthesis of antibiotics', 'Biosynthesis of amino acids',
]

PHARMA_SCORES = {
    'Alkaloids and derivatives': 1.0, 'Alkaloids': 1.0, 'Flavonoids': 0.9,
    'Terpenoids': 0.8, 'Lignans, neolignans and related compounds': 0.8,
    'Phenylpropanoids and polyketides': 0.7, 'Benzenoids': 0.6,
    'Organoheterocyclic compounds': 0.5, 'Organic acids and derivatives': 0.4,
    'Nucleosides, nucleotides, and analogues': 0.4, 'Organic oxygen compounds': 0.3,
    'Lipids and lipid-like molecules': 0.2,
}

# ====================== 核心处理函数 ======================
def read_uploaded_file(uploaded_file) -> pd.DataFrame:
    name = uploaded_file.name
    content = uploaded_file.read()
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
    for df in dfs:
        df.columns = df.columns.str.strip()
    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.sort_values('Vip_plsda', ascending=False)
    merged = merged.drop_duplicates(subset=['Metabolite'], keep='first')
    return merged.reset_index(drop=True)

def filter_differential(df, vip_thresh=1.0, p_thresh=0.05):
    mask = (df['Vip_plsda'] > vip_thresh) & (df['P_value'] < p_thresh)
    return df.loc[mask].copy().reset_index(drop=True)

def annotate_metabolites(diff_df, anno_df):
    return diff_df.merge(anno_df, left_on='Metabolite', right_on='metab', how='left')

def is_pharmacological(row) -> tuple:
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
    results = []
    for _, row in df.iterrows():
        is_pharma, reason = is_pharmacological(row)
        if is_pharma:
            row_dict = row.to_dict()
            row_dict['_pharma_reason'] = reason
            results.append(row_dict)
    return pd.DataFrame(results)

def calculate_abundance(pharma_metabolites: list, diff_df: pd.DataFrame) -> pd.DataFrame:
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
    df = pharma_df.merge(abundance_df, on='Metabolite', how='left')
    vip_max = df['Vip_plsda'].max()
    df['_vip_norm'] = df['Vip_plsda'] / vip_max if vip_max > 0 else 0
    species_cols = [c for c in df.columns if c.endswith('_mean_abundance')]
    df['_max_abund'] = df[species_cols].replace(0, np.nan).max(axis=1)
    df['_min_abund'] = df[species_cols].replace(0, np.nan).min(axis=1)
    df['_diff_ratio'] = df['_max_abund'] / (df['_min_abund'] + 1e-9)
    df['_diff_ratio'] = df['_diff_ratio'].replace([np.inf, -np.inf], np.nan)
    diff_max = df['_diff_ratio'].quantile(0.95)
    df['_diff_norm'] = np.clip(df['_diff_ratio'] / diff_max, 0, 1)
    df['_pharma_score'] = df['super_class'].map(lambda x: PHARMA_SCORES.get(str(x).strip(), 0.3))
    df['_star_score'] = df['_vip_norm'] * 0.3 + df['_diff_norm'].fillna(0) * 0.3 + df['_pharma_score'] * 0.4
    return df.sort_values('_star_score', ascending=False).reset_index(drop=True)

def extract_kegg_ids(df) -> list:
    kegg_ids = []
    for _, row in df.iterrows():
        pid = str(row.get('pathway_id', '-'))
        for p in pid.split(';'):
            p = p.strip()
            if p.startswith('map') or p.startswith('ko'):
                kegg_ids.append(p)
    return list(set(kegg_ids))

def df_to_excel_bytes(df: pd.DataFrame) -> bytes:
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Sheet1')
    output.seek(0)
    return output.getvalue()

# ====================== 可视化函数 ======================
def plot_heatmap(pharma_df: pd.DataFrame, abundance_df: pd.DataFrame, top_n=50):
    from scipy import stats
    species_order = ['LJ', 'BC', 'XY', 'EH', 'MB']
    abundance_cols = [f'{sp}_mean_abundance' for sp in species_order if f'{sp}_mean_abundance' in abundance_df.columns]
    if not abundance_cols:
        return None
    top_metabs = pharma_df.nlargest(top_n, 'Vip_plsda')['Metabolite'].tolist()
    plot_df = abundance_df[abundance_df['Metabolite'].isin(top_metabs)].copy().set_index('Metabolite').loc[top_metabs]
    z_df = plot_df[abundance_cols].apply(lambda x: stats.zscore(x, nan_policy='omit'), axis=0)
    z_df.columns = [c.replace('_mean_abundance', '') for c in z_df.columns]
    z_df.index = [name[:30] + '...' if len(name) > 30 else name for name in z_df.index]
    fig, ax = plt.subplots(figsize=(8, max(10, top_n * 0.25)))
    sns.heatmap(z_df, cmap='RdBu_r', center=0, xticklabels=True, yticklabels=True, cbar_kws={'label': 'Z-score'}, ax=ax)
    ax.set_title(f'Candidate Compounds Heatmap (Top {top_n} by VIP)', fontsize=13)
    ax.set_xlabel('Species'); ax.set_ylabel('Metabolite')
    plt.xticks(rotation=0); plt.yticks(fontsize=7); plt.tight_layout()
    return fig

def plot_boxplot(abundance_df: pd.DataFrame, pharma_df: pd.DataFrame, top_n=20):
    species_order = ['LJ', 'BC', 'XY', 'EH', 'MB']
    abundance_cols = [f'{sp}_mean_abundance' for sp in species_order if f'{sp}_mean_abundance' in abundance_df.columns]
    if not abundance_cols:
        return None
    top_metabs = pharma_df.nlargest(top_n, 'Vip_plsda')['Metabolite'].tolist()
    plot_df = abundance_df[abundance_df['Metabolite'].isin(top_metabs)].copy()
    melt_df = plot_df.melt(id_vars=['Metabolite'], value_vars=abundance_cols, var_name='Species', value_name='Abundance')
    melt_df['Species'] = melt_df['Species'].str.replace('_mean_abundance', '')
    melt_df['Abundance'] = pd.to_numeric(melt_df['Abundance'], errors='coerce')
    melt_df = melt_df.dropna(subset=['Abundance'])
    melt_df['Metabolite_short'] = melt_df['Metabolite'].apply(lambda x: x[:25] + '...' if len(x) > 25 else x)
    fig = px.box(melt_df, x='Species', y='Abundance', color='Species', title=f'Abundance Distribution per Species (Top {top_n} VIP)', points='outliers', hover_data={'Metabolite': True, 'Abundance': ':.4f'})
    fig.update_layout(boxmode='group', showlegend=True, xaxis_title='Species', yaxis_title='Mean Abundance (log2)', font=dict(size=12))
    return fig

def plot_radar(pharma_row: dict, abundance_row: dict) -> go.Figure:
    species_order = ['LJ', 'BC', 'XY', 'EH', 'MB']
    abundance_cols = [f'{sp}_mean_abundance' for sp in species_order if f'{sp}_mean_abundance' in abundance_row]
    vals, labels = [], []
    for col in abundance_cols:
        sp = col.replace('_mean_abundance', '')
        v = abundance_row.get(col)
        if v is not None and not pd.isna(v):
            vals.append(float(v)); labels.append(sp)
        else:
            vals.append(0); labels.append(sp)
    if len(vals) < 3:
        return None
    vals = vals + [vals[0]]; labels = labels + [labels[0]]
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(r=vals, theta=labels, fill='toself', fillcolor='rgba(31,119,180,0.3)', line_color='rgb(31,119,180)', name=pharma_row.get('Metabolite', '')[:20]))
    fig.update_layout(polar=dict(radialaxis=dict(visible=True)), showlegend=False, title=dict(text=pharma_row.get('Metabolite', '')[:40], font=dict(size=12)), height=300)
    return fig

# ====================== Streamlit UI ======================
def main():
    st.markdown("""
    <div style="background: linear-gradient(135deg, #630ed4 0%, #7c3aed 50%, #9f67f5 100%);
                padding: 24px 32px; border-radius: 1.5rem; margin-bottom: 24px;
                box-shadow: 0 8px 32px rgba(99,14,212,0.3);">
        <div style="display: flex; align-items: center; justify-content: space-between;">
            <div style="display: flex; align-items: center; gap: 16px;">
                <div style="font-size: 36px; filter: drop-shadow(0 2px 8px rgba(0,0,0,0.2));">🧬</div>
                <div>
                    <h1 style="color: #ffffff; margin: 0; font-size: 22px; font-weight: 700; letter-spacing: -0.03em; text-shadow: 0 1px 4px rgba(0,0,0,0.15);">差异代谢物药用筛选平台</h1>
                    <p style="color: rgba(255,255,255,0.82); margin: 6px 0 0 0; font-size: 13px; letter-spacing: 0.01em;">Multivariate Stats · Volcano Plot · KEGG Enrichment · Network Pharmacology</p>
                </div>
            </div>
            <div style="background: rgba(255,255,255,0.15); border-radius: 1rem; padding: 8px 16px; backdrop-filter: blur(8px);">
                <span style="color: rgba(255,255,255,0.9); font-size: 12px; font-weight: 600; letter-spacing: 0.05em; text-transform: uppercase;">MetaboLab</span>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # ===== 侧边栏 =====
    st.sidebar.title("Config")
    st.sidebar.markdown("---")
    uploaded_diff = st.sidebar.file_uploader("上传 diff.exp.xls 文件（支持多选）", type=['xls', 'xlsx', 'tsv', 'txt'], accept_multiple_files=True, help="同时选中所有比较组的文件")
    uploaded_anno = st.sidebar.file_uploader("上传 anno.xls 注释文件", type=['xls', 'xlsx', 'tsv', 'txt'], help="包含化合物注释信息的anno.xls文件")
    st.sidebar.markdown("---")
    st.sidebar.subheader("Parameters")
    vip_thresh = st.sidebar.slider("VIP Threshold (Vip_plsda > ?)", 0.5, 3.0, 1.0, 0.1, help="只保留 VIP > 此值的差异代谢物")
    p_thresh = st.sidebar.slider("P-value Threshold (P_value < ?)", 0.01, 0.1, 0.05, 0.005, help="只保留 P_value < 此值的差异代谢物")
    run_analysis = st.sidebar.button("Start Analysis", type="primary", width="stretch")
    st.sidebar.markdown("---")
    st.sidebar.markdown("**Input format**: `.xls`/`.xlsx`/`.tsv` 文件（自动识别）")

    st.title("差异代谢物药用筛选 Web")
    st.markdown("上传差异分析文件 + 注释文件，自动完成筛选、可视化、明星分子精选、KEGG 富集分析。")

    if 'analysis_done' not in st.session_state:
        st.session_state['analysis_done'] = False
        st.session_state['diff_df'] = None
        st.session_state['pharma_df'] = None
        st.session_state['abundance_df'] = None
        st.session_state['annotated_df'] = None
        st.session_state['star_df'] = None

    if run_analysis:
        if not uploaded_diff:
            st.warning("请至少上传一个 diff.exp.xls 文件！")
            return
        if not uploaded_anno:
            st.warning("请上传 anno.xls 注释文件！")
            return

        with st.spinner("Reading files..."):
            diff_dfs = []
            for f in uploaded_diff:
                try:
                    df = read_uploaded_file(f)
                    diff_dfs.append(df)
                    st.sidebar.success(f"Loaded: {f.name} ({len(df)} rows)")
                except Exception as e:
                    st.sidebar.error(f"Failed to read {f.name}: {e}")
                    return

            anno_df = read_uploaded_file(uploaded_anno)

        with st.spinner("Merging diff files..."):
            merged = process_diff_files(diff_dfs)
            st.info(f"Merged: {len(merged)} unique metabolites from {len(diff_dfs)} files")

        with st.spinner("Filtering differential metabolites..."):
            diff_filtered = filter_differential(merged, vip_thresh, p_thresh)
            st.info(f"Differential filter (VIP>{vip_thresh}, P<{p_thresh}): {len(diff_filtered)} metabolites")

        with st.spinner("Annotating metabolites..."):
            annotated = annotate_metabolites(diff_filtered, anno_df)
            named_count = annotated['compound_name'].notna().sum() - (annotated['compound_name'] == '-').sum()
            st.info(f"Annotation: {int(named_count)} named compounds matched")

        with st.spinner("Filtering pharmacological candidates..."):
            pharma = filter_pharmacological(annotated)
            st.info(f"Pharmacological candidates: {len(pharma)}")

        with st.spinner("Calculating abundance..."):
            pharma_metabs = pharma['Metabolite'].tolist()
            abundance = calculate_abundance(pharma_metabs, merged)
            st.info(f"Abundance calculated for {len(abundance)} compounds")

        with st.spinner("Computing star scores..."):
            star = compute_star_scores(pharma, abundance)
            st.info(f"Star scores computed")

        st.session_state['diff_df'] = merged
        st.session_state['diff_filtered'] = diff_filtered
        st.session_state['annotated_df'] = annotated
        st.session_state['pharma_df'] = pharma
        st.session_state['abundance_df'] = abundance
        st.session_state['star_df'] = star
        st.session_state['analysis_done'] = True
        st.session_state['uploaded_diff'] = uploaded_diff
        st.session_state['uploaded_anno'] = uploaded_anno
        st.session_state['vip_thresh'] = vip_thresh
        st.session_state['p_thresh'] = p_thresh
        st.success("Analysis complete!")

    if st.session_state['analysis_done']:
        tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
            "📊 Data Overview & Downloads",
            "📈 Visualization",
            "⭐ Star Molecules",
            "🧬 KEGG Enrichment",
            "💊 Pharma DB",
            "🧮 Multivariate Stats",
            "🔥 Two-Group Comparison",
            "📊 Multi-Group Comparison",
            "🧬 Network Pharmacology",
        ])

        with tab1:
            st.subheader("Data Overview")
            col_stat1, col_stat2, col_stat3, col_stat4 = st.columns(4)
            col_stat1.metric("Total Metabolites", len(st.session_state['diff_df']))
            col_stat2.metric("Differential (VIP>P)", len(st.session_state['diff_filtered']))
            named = st.session_state['annotated_df']['compound_name'].notna() & (st.session_state['annotated_df']['compound_name'] != '-')
            col_stat3.metric("Named Compounds", int(named.sum()))
            col_stat4.metric("Pharmacological", len(st.session_state['pharma_df']))

            st.markdown("---")
            st.subheader("Download Results")
            col_dl1, col_dl2, col_dl3 = st.columns(3)

            diff_out_cols = ['Metabolite', 'Vip_plsda', 'Vip_oplsda', 'P_value', 'fdr', 'FC', 'compound_name', 'formula', 'super_class', 'kegg_second_category']
            diff_out_cols = [c for c in diff_out_cols if c in st.session_state['annotated_df'].columns]
            diff_out = st.session_state['annotated_df'][diff_out_cols].drop_duplicates(subset=['Metabolite'])
            bytes1 = df_to_excel_bytes(diff_out)
            col_dl1.download_button(label="Download: Total Differential Metabolites", data=bytes1, file_name="total_differential_metabolites.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")
            col_dl1.caption(f"{len(diff_out)} rows")

            out2_cols = ['Metabolite', 'compound_name', 'formula', 'hmdb_id', 'super_class', 'class', 'sub_class', 'compound_first_category', 'compound_second_category', 'kegg_first_category', 'kegg_second_category', 'pathway_id', 'Vip_plsda', 'P_value', 'FC', '_pharma_reason']
            out2_cols = [c for c in out2_cols if c in st.session_state['pharma_df'].columns]
            pharma_out = st.session_state['pharma_df'][out2_cols].copy().rename(columns={'_pharma_reason': '筛选依据'})
            bytes2 = df_to_excel_bytes(pharma_out)
            col_dl2.download_button(label="Download: Pharmacological Candidates", data=bytes2, file_name="pharmacological_candidates.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")
            col_dl2.caption(f"{len(pharma_out)} rows")

            bytes3 = df_to_excel_bytes(st.session_state['abundance_df'])
            col_dl3.download_button(label="Download: Five-Species Abundance", data=bytes3, file_name="five_species_abundance.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")
            col_dl3.caption(f"{len(st.session_state['abundance_df'])} rows")

            st.markdown("---")
            if generate_word_report is not None:
                st.subheader("Download Analysis Report")
                if st.button("Generate Word Report", type="secondary", width="stretch"):
                    with st.spinner("Generating Word report..."):
                        try:
                            from io import BytesIO
                            annotated_df = st.session_state['annotated_df']
                            pharma_df = st.session_state['pharma_df']
                            star_df = st.session_state.get('star_df', pharma_df)
                            abundance_df = st.session_state['abundance_df']
                            params = {
                                'vip_thresh': st.session_state.get('vip_thresh', 1.0),
                                'p_thresh': st.session_state.get('p_thresh', 0.05),
                                'n_diff': len(st.session_state['diff_filtered']),
                                'n_named': int((annotated_df['compound_name'].notna() & (annotated_df['compound_name'] != '-')).sum()),
                                'n_pharma': len(pharma_df),
                                'n_star': len(star_df),
                            }
                            doc_bytes = generate_word_report(
                                annotated_df=annotated_df, pharma_df=pharma_df, star_df=star_df,
                                abundance_df=abundance_df, params=params,
                                fig_heatmap_bytes=None, fig_boxplot_bytes=None,
                                fig_vip_bytes=None, fig_bubble_bytes=None, fig_enr_bar_bytes=None,
                                enr_df=st.session_state.get('enr_df'),
                            )
                        except Exception as e:
                            doc_bytes = None
                            st.error(f"Report generation failed: {e}")
                        if doc_bytes:
                            st.success("Report generated successfully!")
                            st.download_button("Download Word Report (.docx)", data=doc_bytes, file_name=f"metabolite_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M')}.docx", mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document", width="stretch")

            st.markdown("---")
            st.subheader("Preview: Pharmacological Candidates")
            display_cols = ['Metabolite', 'compound_name', 'super_class', 'Vip_plsda', 'P_value', 'FC', '_pharma_reason']
            display_cols = [c for c in display_cols if c in st.session_state['pharma_df'].columns]
            st.dataframe(st.session_state['pharma_df'][display_cols].head(30), width="stretch", height=400)

        with tab2:
            st.subheader("Visualization")
            viz_col1, viz_col2 = st.columns([1, 1])
            with viz_col1:
                st.markdown("**Clustered Heatmap (Top N by VIP)**")
                hm_top_n = st.slider("Select Top N compounds for heatmap:", 10, 100, 30, 5, key="hm_top")
                fig_hm = plot_heatmap(st.session_state['pharma_df'], st.session_state['abundance_df'], top_n=hm_top_n)
                if fig_hm:
                    st.pyplot(fig_hm)
                    buf_fig = io.BytesIO()
                    fig_hm.savefig(buf_fig, format='png', dpi=150, bbox_inches='tight')
                    buf_fig.seek(0)
                    st.download_button("Download Heatmap PNG", data=buf_fig, file_name="heatmap.png", mime="image/png", width="stretch")
                else:
                    st.warning("Insufficient data for heatmap.")
            with viz_col2:
                st.markdown("**Boxplot by Species (Top N by VIP)**")
                bx_top_n = st.slider("Select Top N compounds for boxplot:", 5, 50, 15, 5, key="bx_top")
                fig_bx = plot_boxplot(st.session_state['abundance_df'], st.session_state['pharma_df'], top_n=bx_top_n)
                if fig_bx:
                    st.plotly_chart(fig_bx, width="stretch")
                else:
                    st.warning("Insufficient data for boxplot.")
            st.markdown("---")
            st.subheader("VIP Score Distribution")
            fig_vip = px.histogram(st.session_state['annotated_df'], x='Vip_plsda', nbins=50, title='VIP Score Distribution (Differential Metabolites)', labels={'Vip_plsda': 'VIP (PLS-DA)', 'count': 'Count'}, color_discrete_sequence=['#4C78A8'])
            fig_vip.add_vline(x=vip_thresh, line_dash='dash', line_color='red', annotation_text=f'VIP threshold={vip_thresh}')
            fig_vip.update_layout(height=400)
            st.plotly_chart(fig_vip, width="stretch")

        with tab3:
            st.subheader("Star Molecule Selection")
            if 'star_df_v2' in st.session_state:
                star_df = st.session_state['star_df_v2'].copy()
                score_col = '_star_score_v2'
                st.info("已更新为含药理证据的综合评分 (v2)")
            else:
                star_df = st.session_state['star_df'].copy()
                score_col = '_star_score'
            all_super_classes = sorted(star_df['super_class'].dropna().unique().tolist())
            col_f1, col_f2, col_f3 = st.columns([1, 1, 1])
            with col_f1:
                selected_classes = st.multiselect("Pharmacological Category", options=all_super_classes, default=['Alkaloids and derivatives', 'Phenylpropanoids and polyketides'], format_func=lambda x: x if x else 'All')
            with col_f2:
                min_vip = st.slider("Minimum VIP", 0.5, 5.0, 1.0, 0.1)
            with col_f3:
                min_diff_ratio = st.slider("Minimum Species Diff Ratio (max/min)", 1.0, 100.0, 1.0, 0.5)
            filtered_star = star_df[star_df['Vip_plsda'] >= min_vip].copy()
            if selected_classes:
                filtered_star = filtered_star[filtered_star['super_class'].isin(selected_classes)]
            if '_diff_ratio' in filtered_star.columns:
                filtered_star = filtered_star[filtered_star['_diff_ratio'] >= min_diff_ratio]
            st.markdown(f"Filtered: **{len(filtered_star)}** / {len(star_df)} candidates")
            st.markdown("---")
            st.subheader("Top 10 Star Molecules (by Composite Score)")
            top10 = filtered_star.head(10)
            for idx, row in top10.iterrows():
                rank = list(top10['Metabolite']).index(row['Metabolite']) + 1
                with st.container():
                    c1, c2 = st.columns([3, 1])
                    with c1:
                        metab_name = row.get('Metabolite', '')
                        compound = row.get('compound_name', '-')
                        sc = row.get('super_class', '-')
                        vip = row.get('Vip_plsda', 0)
                        diff_r = row.get('_diff_ratio', 0)
                        score = row.get(score_col, row.get('_star_score', 0))
                        reason = row.get('_pharma_reason', '')
                        st.markdown(f"**#{rank} {metab_name}** | VIP={vip:.2f} | Diff={diff_r:.1f}x | Score={score:.3f} | *{sc}* | {reason}")
                    with c2:
                        search_name = str(compound) if compound and compound != '-' else str(metab_name)
                        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={search_name.replace(' ', '+')}"
                        scholar_url = f"https://scholar.google.com/scholar?q={search_name.replace(' ', '+')}"
                        st.markdown(f"[PubMed]({pubmed_url}) · [Google Scholar]({scholar_url})", unsafe_allow_html=True)
                    st.divider()
            st.markdown("---")
            st.subheader("Full Filtered List (with Scores)")
            table_cols = ['Metabolite', 'compound_name', 'super_class', 'Vip_plsda', '_diff_ratio', score_col, '_vip_norm', '_diff_norm', '_pharma_score', 'LJ_mean_abundance', 'BC_mean_abundance', 'XY_mean_abundance', 'EH_mean_abundance', 'MB_mean_abundance']
            table_cols = [c for c in table_cols if c in filtered_star.columns]
            st.dataframe(filtered_star[table_cols].style.format({'Vip_plsda': '{:.4f}', '_diff_ratio': '{:.2f}', score_col: '{:.3f}', '_vip_norm': '{:.3f}', '_diff_norm': '{:.3f}', '_pharma_score': '{:.2f}', 'LJ_mean_abundance': '{:.4f}', 'BC_mean_abundance': '{:.4f}', 'XY_mean_abundance': '{:.4f}', 'EH_mean_abundance': '{:.4f}', 'MB_mean_abundance': '{:.4f}'}, na_rep='-'), width="stretch", height=500)
            bytes_star = df_to_excel_bytes(filtered_star[table_cols])
            st.download_button("Download Star Molecules Table", data=bytes_star, file_name="star_molecules.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")

        with tab4:
            st.subheader("KEGG Pathway Enrichment Analysis")
            st.markdown("基于注释数据计算 KEGG 通路过representation分析（超几何检验），无需联网。")
            if st.button("Run Pathway Enrichment", type="primary", width="stretch"):
                pharma_df = st.session_state['pharma_df']
                annotated_df = st.session_state['annotated_df']
                metab_to_pathways = {}
                for _, row in annotated_df.iterrows():
                    metab = row['Metabolite']
                    pid = str(row.get('pathway_id', '-'))
                    if pid not in ('-', 'nan', '') and pd.notna(row.get('pathway_id')):
                        pathways = [p.strip() for p in pid.split(';') if p.strip().startswith('map')]
                        if pathways:
                            metab_to_pathways[metab] = pathways
                bg_mets = set(metab_to_pathways.keys())
                N_bg = len(bg_mets)
                pathway_bg_count = {}
                for pathways in metab_to_pathways.values():
                    for p in pathways:
                        pathway_bg_count[p] = pathway_bg_count.get(p, 0) + 1
                pharma_mets = set(pharma_df['Metabolite'].unique())
                n_cand = len(pharma_mets & set(metab_to_pathways.keys()))
                st.info(f"背景代谢物（注释到KEGG通路）: {N_bg} | 候选化合物（注释到KEGG通路）: {n_cand}")
                enr_results = []
                for pathway_id, K in pathway_bg_count.items():
                    hit_mets = [m for m in pharma_mets if m in metab_to_pathways and pathway_id in metab_to_pathways[m]]
                    k = len(hit_mets)
                    if k > 0:
                        pval = stats.hypergeom.sf(k - 1, N_bg, K, n_cand)
                        hit_names = '; '.join(sorted(hit_mets)[:8])
                        if len(hit_mets) > 8:
                            hit_names += f' ... (+{len(hit_mets) - 8} more)'
                        enr_results.append({'KEGG_ID': pathway_id, 'KEGG_URL': f'https://www.genome.jp/kegg/mapper/?org_name=hsa&mode=compound&pathway={pathway_id}', 'Pathway_Name': _PATHWAY_NAME_MAP.get(pathway_id, pathway_id), 'Candidate_Hits': k, 'Total_in_Pathway': K, 'Candidate_Total': n_cand, 'Rich_Factor': round(k / K, 3) if K > 0 else 0, 'P_value': pval, 'Hit_Metabolites': hit_names})
                enr_df = pd.DataFrame(enr_results)
                if len(enr_df) > 0:
                    enr_df = enr_df.sort_values('P_value').reset_index(drop=True)
                    n_tests = len(enr_df)
                    adj_pvals = []
                    for i, (_, row) in enumerate(enr_df.iterrows()):
                        rank = i + 1
                        adj_p = min(row['P_value'] * n_tests / rank, 1.0)
                        adj_pvals.append(adj_p)
                    enr_df['P_value_BH'] = adj_pvals
                    enr_df['-log10P'] = -np.log10(enr_df['P_value'].replace(0, 1e-10))
                    sig_df = enr_df[enr_df['P_value'] < 0.05].copy()
                    st.success(f"显著富集通路: {len(sig_df)} 条 (P < 0.05)")
                    st.session_state['enr_df'] = sig_df
                    if len(sig_df) > 0:
                        st.markdown("**Pathway Enrichment Bubble Chart (P < 0.05)**")
                        plot_df = sig_df.head(30).copy()
                        plot_df['label'] = plot_df['KEGG_ID'] + ': ' + plot_df['Pathway_Name']
                        plot_df['label'] = plot_df['label'].apply(lambda x: x[:55] + '...' if len(str(x)) > 55 else x)
                        fig_bubble = px.scatter(plot_df, x='Rich_Factor', y='label', size='Candidate_Hits', color='-log10P', size_max=40, color_continuous_scale='RdYlBu_r', title='KEGG Pathway Enrichment Bubble Chart (Top 30 by P-value)', labels={'Rich_Factor': 'Rich Factor', 'label': 'Pathway (KEGG ID: Name)', 'Candidate_Hits': '# Compounds', '-log10P': '-log10(P-value)'}, hover_data={'KEGG_ID': True, 'Candidate_Hits': True, 'P_value': ':.2e'})
                        fig_bubble.update_layout(height=max(400, len(plot_df) * 18))
                        fig_bubble.update_yaxes(tickfont=dict(size=9))
                        st.plotly_chart(fig_bubble, width="stretch")
                        st.markdown("**Top Pathways by Hit Count**")
                        bar_df = sig_df.nlargest(20, 'Candidate_Hits')
                        bar_df['label'] = bar_df['KEGG_ID'] + ': ' + bar_df['Pathway_Name']
                        bar_df['label'] = bar_df['label'].apply(lambda x: x[:55] + '...' if len(str(x)) > 55 else x)
                        fig_bar = px.bar(bar_df, x='Candidate_Hits', y='label', orientation='h', color='-log10P', color_continuous_scale='Viridis', title='Top 20 Enriched Pathways by Hit Count', labels={'Candidate_Hits': '# Compounds in Pathway', 'label': 'Pathway', '-log10P': '-log10(P-value)'}, hover_data={'Rich_Factor': True, 'P_value': ':.2e'})
                        fig_bar.update_layout(height=max(400, len(bar_df) * 18), yaxis={'autorange': 'reversed'})
                        fig_bar.update_yaxes(tickfont=dict(size=9))
                        st.plotly_chart(fig_bar, width="stretch")
                        display_cols = ['KEGG_ID', 'Pathway_Name', 'Candidate_Hits', 'Rich_Factor', 'P_value', 'P_value_BH', 'Hit_Metabolites']
                        st.dataframe(sig_df[display_cols].style.format({'Rich_Factor': '{:.3f}', 'P_value': '{:.2e}', 'P_value_BH': '{:.2e}'}, na_rep='-'), width="stretch", height=400)
                        bytes_enr = df_to_excel_bytes(sig_df[display_cols])
                        st.download_button("Download Enrichment Results (P<0.05)", data=bytes_enr, file_name="pathway_enrichment_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")
                else:
                    st.warning("候选化合物中未找到带 KEGG pathway_id 的注释。")
            st.markdown("---")
            st.subheader("Pathway Coverage Overview")
            if st.session_state.get('pharma_df') is not None:
                pharma_df = st.session_state['pharma_df']
                all_pathways = []
                for pid in pharma_df['pathway_id'].dropna():
                    if str(pid) not in ('-', 'nan', ''):
                        for p in str(pid).split(';'):
                            p = p.strip()
                            if p.startswith('map') or p.startswith('ko'):
                                all_pathways.append(p)
                if all_pathways:
                    pw_counts = pd.Series(all_pathways).value_counts()
                    pw_display = pw_counts.head(30).reset_index()
                    pw_display.columns = ['KEGG_ID', 'Count']
                    pw_display['Pathway_Name'] = pw_display['KEGG_ID'].map(_PATHWAY_NAME_MAP).fillna(pw_display['KEGG_ID'])
                    fig_dist = px.bar(pw_display, x='KEGG_ID', y='Count', color='Count', color_continuous_scale='Blues', title='KEGG Pathway ID Distribution (Top 30)', labels={'KEGG_ID': 'KEGG Pathway ID', 'Count': 'Compound Count'}, text='Pathway_Name')
                    fig_dist.update_layout(height=450, xaxis={'tickangle': 45})
                    st.plotly_chart(fig_dist, width="stretch")

        with tab5:
            st.subheader("Pharmacology Database Integration (Auto-Download)")
            if match_pharma_online is None:
                st.error("pharma_cache.py 模块未正确加载。请确保 pharma_cache.py 与 app.py 在同一目录下。")
            else:
                # 显示本地缓存状态提示
                cache_info = get_cache_status()
                if cache_info['disk_cache_valid']:
                    st.success(f"✅ 检测到本地 DrugCentral 缓存（共 {cache_info.get('disk_cache_targets_count', 0)} 条记录），将直接使用，无需重新下载")
                else:
                    st.info("ℹ️ 首次使用将自动下载 DrugCentral 数据库（约 1-2 分钟），之后自动使用本地缓存")

                st.markdown("#### 自动在线查询（DrugCentral + PubChem）")
                online_query_btn = st.button("Start Online Pharma Query (Auto-Download)", type="primary", width="stretch", help="首次运行时会自动下载 DrugCentral 数据库（约 1-2 分钟），之后直接使用缓存")
                if online_query_btn:
                    # 根据缓存状态动态显示 spinner 提示
                    if cache_info['disk_cache_valid']:
                        spinner_msg = "正在从本地缓存加载 DrugCentral 数据..."
                    else:
                        spinner_msg = "首次使用，正在下载 DrugCentral 数据库（仅下载一次）..."
                    with st.spinner(spinner_msg):
                        metabolites = st.session_state['pharma_df']['Metabolite'].tolist()
                        progress_bar = st.progress(0)
                        def progress_callback(current, total):
                            progress_bar.progress(int(current / total * 100))
                        pharma_match_df = match_pharma_online(metabolites=metabolites, progress_callback=progress_callback)
                        progress_bar.empty()
                        st.session_state['pharma_match_df'] = pharma_match_df
                        n_dc = pharma_match_df['DrugCentral_Targets'].apply(lambda x: x != '-' and len(str(x)) > 0 if pd.notna(x) else False).sum()
                        n_pc = (pharma_match_df['PubChem_BioActivity_Count'] > 0).sum()
                        n_active = (pharma_match_df['Pharma_Evidence_Score'] > 0).sum()
                        st.success(f"查询完成！DrugCentral命中: {n_dc} | PubChem生物活性: {n_pc} | 有药理证据: {n_active}/{len(pharma_match_df)}")
                        if 'star_df' in st.session_state and len(st.session_state['star_df']) > 0:
                            with st.spinner("Recalculating star scores with pharma evidence..."):
                                star_df = st.session_state['star_df'].copy()
                                score_map = pharma_match_df.set_index('Metabolite')['Pharma_Evidence_Score'].to_dict()
                                star_df['_pharma_evidence'] = star_df['Metabolite'].map(score_map).fillna(0)
                                pe_max = star_df['_pharma_evidence'].max()
                                star_df['_pharma_evidence_norm'] = star_df['_pharma_evidence'] / pe_max if pe_max > 0 else 0
                                star_df['_star_score_v2'] = star_df['_vip_norm'].fillna(0) * 0.25 + star_df['_diff_norm'].fillna(0) * 0.25 + star_df['_pharma_score'].fillna(0) * 0.30 + star_df['_pharma_evidence_norm'].fillna(0) * 0.20
                                star_df = star_df.sort_values('_star_score_v2', ascending=False).reset_index(drop=True)
                                st.session_state['star_df_v2'] = star_df
                                st.success("明星分子评分已更新（含药理证据维度）！请切换到 Star Molecules 标签页查看。")

                # ---- 可选：本地文件上传区（作为补充）----
                st.markdown("---")
                st.markdown("#### Optional: Local DB Files / Directory (TCMSP / DrugBank / TTD)")
                if match_pharma_db is not None:
                    tcmsp_dir = st.text_input(
                        "TCMSP 数据目录路径",
                        value="",
                        placeholder=r"C:\Users\tyf\.qclaw\workspace\tcmsp_data",
                        help="填入 TCMSP 数据文件所在目录，程序将自动扫描并加载所有 CSV 文件"
                    ).strip()
                    tcmsp_file = st.file_uploader("TCMSP 单文件", type=['xlsx', 'xls', 'csv'], help="与目录路径二选一")
                    db_col1, db_col2 = st.columns(2)
                    with db_col1:
                        drugbank_file = st.file_uploader("DrugBank File", type=['csv', 'tsv', 'xml'], help="DrugBank 数据库")
                    with db_col2:
                        ttd_file = st.file_uploader("TTD File", type=['xlsx', 'xls', 'csv', 'tsv'], help="TTD (Therapeutic Target Database)")
                    local_query_btn = st.button("Run Local DB Query", width="stretch")
                    if local_query_btn:
                        with st.spinner("Loading local database files..."):
                            if tcmsp_dir:
                                try:
                                    tcmsp_data = load_tcmsp_from_dir(tcmsp_dir)
                                    st.info(f"TCMSP 目录加载: {len(tcmsp_data)} 化合物")
                                except Exception as e:
                                    st.warning(f"TCMSP 目录加载失败: {e}，尝试单文件方式")
                                    tcmsp_data = load_tcmsp(tcmsp_file) if tcmsp_file else {}
                            else:
                                tcmsp_data = load_tcmsp(tcmsp_file) if tcmsp_file else {}
                            drugbank_data = load_drugbank(drugbank_file) if drugbank_file else {}
                            ttd_data = load_ttd(ttd_file) if ttd_file else {}
                            st.info(f"TCMSP: {len(tcmsp_data)} | DrugBank: {len(drugbank_data)} | TTD: {len(ttd_data)}")
                        with st.spinner("Matching compounds..."):
                            metabolites = st.session_state['pharma_df']['Metabolite'].tolist()
                            progress_bar2 = st.progress(0)
                            def pbar(current, total):
                                progress_bar2.progress(int(current / total * 100))
                            pharma_match_df = match_pharma_db(
                                metabolites=metabolites,
                                tcmsp_data=tcmsp_data,
                                drugbank_data=drugbank_data,
                                ttd_data=ttd_data,
                                progress_callback=pbar,
                            )
                            progress_bar2.empty()
                            st.session_state['pharma_match_df'] = pharma_match_df
                            n_tcmsp = pharma_match_df['TCMSP_OB'].notna().sum()
                            n_db = (pharma_match_df['DrugBank_Targets'] != '-').sum()
                            n_ttd = (pharma_match_df['TTD_Targets'] != '-').sum()
                            n_active = (pharma_match_df['Pharma_Evidence_Score'] > 0).sum()
                            st.success(f"本地匹配完成！TCMSP: {n_tcmsp} | DrugBank: {n_db} | TTD: {n_ttd} | 有证据: {n_active}/{len(pharma_match_df)}")
                            if 'star_df' in st.session_state and len(st.session_state['star_df']) > 0:
                                with st.spinner("Recalculating star scores..."):
                                    star_df = st.session_state['star_df'].copy()
                                    score_map = pharma_match_df.set_index('Metabolite')['Pharma_Evidence_Score'].to_dict()
                                    star_df['_pharma_evidence'] = star_df['Metabolite'].map(score_map).fillna(0)
                                    pe_max = star_df['_pharma_evidence'].max()
                                    star_df['_pharma_evidence_norm'] = star_df['_pharma_evidence'] / pe_max if pe_max > 0 else 0
                                    star_df['_star_score_v2'] = star_df['_vip_norm'].fillna(0) * 0.25 + star_df['_diff_norm'].fillna(0) * 0.25 + star_df['_pharma_score'].fillna(0) * 0.30 + star_df['_pharma_evidence_norm'].fillna(0) * 0.20
                                    star_df = star_df.sort_values('_star_score_v2', ascending=False).reset_index(drop=True)
                                    st.session_state['star_df_v2'] = star_df
                                    st.success("明星分子评分已更新！请切换到 Star Molecules 标签页查看。")

                if 'pharma_match_df' in st.session_state and len(st.session_state['pharma_match_df']) > 0:
                    pm_df = st.session_state['pharma_match_df']
                    is_online = 'Data_Sources' in pm_df.columns
                    st.markdown("---")
                    st.subheader("Pharmacology Match Results")
                    pc1, pc2, pc3, pc4 = st.columns(4)
                    if is_online:
                        n_dc = pm_df['DrugCentral_Targets'].apply(lambda x: x != '-' and len(str(x)) > 0 if pd.notna(x) else False).sum()
                        n_pc = (pm_df['PubChem_BioActivity_Count'] > 0).sum()
                        n_active = (pm_df['Pharma_Evidence_Score'] > 0).sum()
                        pc1.metric("DrugCentral Hits", f"{n_dc}")
                        pc2.metric("PubChem BioActive", f"{n_pc}")
                        pc3.metric("PubChem CID", f"{pm_df['PubChem_CID'].notna().sum()}")
                        pc4.metric("With Evidence", f"{n_active}/{len(pm_df)}")
                    else:
                        n_tcmsp = pm_df['TCMSP_OB'].notna().sum()
                        n_db = (pm_df['DrugBank_Targets'] != '-').sum()
                        n_ttd = (pm_df['TTD_Targets'] != '-').sum()
                        n_active = (pm_df['Pharma_Evidence_Score'] > 0).sum()
                        pc1.metric("TCMSP Hits", f"{n_tcmsp}")
                        pc2.metric("DrugBank Hits", f"{n_db}")
                        pc3.metric("TTD Hits", f"{n_ttd}")
                        pc4.metric("With Evidence", f"{n_active}/{len(pm_df)}")
                    st.markdown("**Pharma Evidence Score Distribution**")
                    fig_ev = px.histogram(pm_df, x='Pharma_Evidence_Score', nbins=20, title='Pharmacological Evidence Score Distribution', labels={'Pharma_Evidence_Score': 'Evidence Score', 'count': 'Count'}, color_discrete_sequence=['#2ECC71'])
                    fig_ev.update_layout(height=350)
                    st.plotly_chart(fig_ev, width="stretch")
                    if is_online:
                        display_cols = ['Metabolite', 'Data_Sources', 'DrugCentral_Targets', 'DrugCentral_Diseases', 'PubChem_CID', 'PubChem_XLogP', 'PubChem_BioActivity_Count', 'PubChem_Targets', 'PubChem_Activities', 'Pharma_Evidence_Score', 'Pharma_Evidence']
                    else:
                        display_cols = ['Metabolite', 'TCMSP_OB', 'TCMSP_DL', 'TCMSP_Targets', 'TCMSP_Diseases', 'DrugBank_Targets', 'DrugBank_Indications', 'TTD_Targets', 'TTD_Diseases', 'Pharma_Evidence_Score', 'Pharma_Evidence']
                    avail_cols = [c for c in display_cols if c in pm_df.columns]
                    st.dataframe(pm_df[avail_cols].style.format({'PubChem_XLogP': lambda x: f'{x:.2f}' if pd.notna(x) else '-', 'PubChem_BioActivity_Count': lambda x: int(x) if pd.notna(x) else 0, 'Pharma_Evidence_Score': '{:.2f}', 'TCMSP_OB': lambda x: f'{x:.1f}' if pd.notna(x) else '-', 'TCMSP_DL': lambda x: f'{x:.3f}' if pd.notna(x) else '-'}, na_rep='-'), width="stretch", height=500)
                    match_bytes = io.BytesIO()
                    with pd.ExcelWriter(match_bytes, engine='openpyxl') as writer:
                        pm_df.to_excel(writer, index=False, sheet_name='Pharma_Match')
                    match_bytes.seek(0)
                    st.download_button("Download Pharma Match Results", data=match_bytes.getvalue(), file_name="pharma_database_match.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")

        with tab6:
            st.subheader("Multivariate Statistical Analysis (PCA / PLS-DA / OPLS-DA)")
            if detect_species_columns_from_diff is None:
                st.error("multivariate_stats 模块未正确加载，请检查是否安装scikit-learn")
            else:
                diff_df = st.session_state.get('diff_df')
                if diff_df is None:
                    st.warning("请先在侧边栏上传文件并点击 Start Analysis")
                else:
                    species_cols = detect_species_columns(diff_df)
                    col_ctrl1, col_ctrl2, col_ctrl3 = st.columns([1, 1, 1])
                    with col_ctrl1:
                        analysis_type = st.selectbox("Analysis Type", ["PCA", "PLS-DA", "OPLS-DA"])
                    with col_ctrl2:
                        group_options = list(species_cols.keys())
                        selected_groups = st.multiselect("Select Groups", group_options, default=group_options[:2] if len(group_options) >= 2 else group_options)
                    with col_ctrl3:
                        scale_method = st.selectbox("Data Scaling", ["UV", "Par"])
                    col_ctrl4, col_ctrl5 = st.columns([1, 1])
                    with col_ctrl4:
                        confidence = st.slider("Confidence Level", 0.80, 0.99, 0.95, 0.01, key="conf_lvl")
                    with col_ctrl5:
                        n_perms = st.number_input("Permutation Tests", 50, 500, 200, 50)
                    run_multivar_btn = st.button("Run Multivariate Analysis", type="primary")
                    if run_multivar_btn:
                        if len(selected_groups) < 2:
                            st.warning("请至少选择2个组进行分析")
                        else:
                            selected_cols = []
                            for sp in selected_groups:
                                selected_cols.extend(species_cols.get(sp, []))
                            available_cols = [c for c in selected_cols if c in diff_df.columns]
                            if len(available_cols) < 3:
                                st.error("有效样本数不足")
                            else:
                                mat_data = diff_df[['Metabolite'] + available_cols].copy().set_index('Metabolite').T
                                mat_data.index.name = 'Sample'
                                groups = []
                                for idx in mat_data.index:
                                    for sp, cols in species_cols.items():
                                        if idx in cols:
                                            groups.append(sp)
                                            break
                                    else:
                                        groups.append('Unknown')
                                mat_data['Group'] = groups
                                st.session_state['_mat_data_groups'] = dict(zip(mat_data.index, mat_data['Group']))
                                with st.spinner("Running analysis..."):
                                    try:
                                        X = mat_data.drop(columns=['Group']).values
                                        if scale_method == 'UV':
                                            X = X / np.std(X, axis=0, ddof=1, keepdims=True)
                                        elif scale_method == 'Par':
                                            sd = np.std(X, axis=0, ddof=1, keepdims=True)
                                            sqrt_n = np.sqrt(np.sum(X**2, axis=0) / X.shape[0])
                                            X = X / (sd / (sqrt_n + 1e-9) + 1e-9)
                                        X_df = pd.DataFrame(X, index=mat_data.index, columns=mat_data.drop(columns=['Group']).columns)
                                        group_labels = mat_data['Group'].values
                                        if analysis_type == "PCA":
                                            scores_df, loadings_df, var_explained = run_pca(X_df, n_components=2)
                                            st.session_state['pca_results'] = {'scores': scores_df, 'loadings': loadings_df, 'variance': var_explained, 'groups': species_cols, 'confidence': confidence}
                                            st.success("PCA analysis completed!")
                                        elif analysis_type == "PLS-DA":
                                            scores_df, loadings_df, vip_dict, var_explained = run_plsda(X_df, group_labels, n_components=2)
                                            st.session_state['plsda_results'] = {'scores': scores_df, 'loadings': loadings_df, 'vip': vip_dict, 'variance': var_explained, 'groups': species_cols, 'confidence': confidence}
                                            st.success("PLS-DA analysis completed!")
                                        elif analysis_type == "OPLS-DA":
                                            scores_df, loadings_df, vip_dict = run_oplsda(X_df, group_labels, n_components=2)
                                            st.session_state['oplsda_results'] = {'scores': scores_df, 'loadings': loadings_df, 'vip': vip_dict, 'groups': species_cols, 'confidence': confidence}
                                            st.success("OPLS-DA analysis completed!")
                                    except Exception as e:
                                        st.error(f"Analysis failed: {e}")
                    _group_map = st.session_state.get('_mat_data_groups', {})
                    st.markdown("---")
                    if analysis_type == "PCA" and 'pca_results' in st.session_state:
                        res = st.session_state['pca_results']
                        scores_df = res['scores'].copy()
                        scores_df['Group'] = [_group_map.get(s, 'Unknown') for s in scores_df['Sample']]
                        col_fig1, col_fig2 = st.columns([1, 1])
                        with col_fig1:
                            fig_scores = plot_scores_scatter(scores_df, scores_df['Group'], f'{analysis_type} Scores Plot', res['variance'][0], res['variance'][1], GROUP_COLORS, 0.15)
                            st.pyplot(fig_scores)
                            buf = io.BytesIO(); fig_scores.savefig(buf, format='png', dpi=150, bbox_inches='tight'); buf.seek(0)
                            st.download_button("Download Scores Plot", data=buf, file_name="scores_plot.png", mime="image/png")
                        with col_fig2:
                            fig_loadings = plot_loading_scatter(res['loadings'], top_n=50, title=f'{analysis_type} Loading Plot')
                            st.pyplot(fig_loadings)
                            buf = io.BytesIO(); fig_loadings.savefig(buf, format='png', dpi=150, bbox_inches='tight'); buf.seek(0)
                            st.download_button("Download Loading Plot", data=buf, file_name="loading_plot.png", mime="image/png")
                    elif analysis_type == "PLS-DA" and 'plsda_results' in st.session_state:
                        res = st.session_state['plsda_results']
                        scores_df = res['scores'].copy()
                        scores_df['Group'] = [_group_map.get(s, 'Unknown') for s in scores_df['Sample']]
                        col_fig1, col_fig2 = st.columns([1, 1])
                        with col_fig1:
                            fig_scores = plot_scores_scatter(scores_df, scores_df['Group'], f'{analysis_type} Scores Plot', res['variance'][0] if len(res['variance']) > 0 else 0.1, res['variance'][1] if len(res['variance']) > 1 else 0.1, GROUP_COLORS, 0.15)
                            st.pyplot(fig_scores)
                            buf = io.BytesIO(); fig_scores.savefig(buf, format='png', dpi=150, bbox_inches='tight'); buf.seek(0)
                            st.download_button("Download Scores Plot", data=buf, file_name="scores_plot.png", mime="image/png")
                        with col_fig2:
                            fig_loadings = plot_loading_scatter(res['loadings'], top_n=50, title=f'{analysis_type} Loading Plot')
                            st.pyplot(fig_loadings)
                            buf = io.BytesIO(); fig_loadings.savefig(buf, format='png', dpi=150, bbox_inches='tight'); buf.seek(0)
                            st.download_button("Download Loading Plot", data=buf, file_name="loading_plot.png", mime="image/png")
                        st.markdown("**Top 20 Metabolites by VIP (PLS-DA)**")
                        vip_df = pd.DataFrame(list(res['vip'].items()), columns=['Metabolite', 'VIP']).sort_values('VIP', ascending=False).head(20)
                        st.dataframe(vip_df, height=300)
                    elif analysis_type == "OPLS-DA" and 'oplsda_results' in st.session_state:
                        res = st.session_state['oplsda_results']
                        scores_df = res['scores'].copy()
                        scores_df['Group'] = [_group_map.get(s, 'Unknown') for s in scores_df['Sample']]
                        fig_scores = plot_scores_scatter(scores_df, scores_df['Group'], f'{analysis_type} Scores Plot', 0.3, 0.2, GROUP_COLORS, 0.15)
                        st.pyplot(fig_scores)
                        buf = io.BytesIO(); fig_scores.savefig(buf, format='png', dpi=150, bbox_inches='tight'); buf.seek(0)
                        st.download_button("Download Scores Plot", data=buf, file_name="scores_plot.png", mime="image/png")
                        st.markdown("**Top 20 Metabolites by VIP (OPLS-DA)**")
                        vip_df = pd.DataFrame(list(res['vip'].items()), columns=['Metabolite', 'VIP']).sort_values('VIP', ascending=False).head(20)
                        st.dataframe(vip_df, height=300)
                    if analysis_type in ["PLS-DA", "OPLS-DA"]:
                        st.markdown("---")
                        st.subheader("Permutation Test")
                        perm_btn = st.button("Run Permutation Test (200 permutations)", type="secondary")
                        if perm_btn:
                            with st.spinner("Running permutation test..."):
                                try:
                                    p_val, perm_df, r2_orig = permutation_test_plsda(X_df.values, group_labels, n_perms=n_perms)
                                    fig_perm = plot_permutation(perm_df, f'Permutation Test ({n_perms} permutations)')
                                    st.pyplot(fig_perm)
                                    st.info(f"Permutation p-value: {p_val:.4f} | Original R²Y: {r2_orig:.4f}")
                                    buf = io.BytesIO(); fig_perm.savefig(buf, format='png', dpi=150, bbox_inches='tight'); buf.seek(0)
                                    st.download_button("Download Permutation Plot", data=buf, file_name="permutation_test.png", mime="image/png")
                                except Exception as e:
                                    st.error(f"Permutation test failed: {e}")

        with tab7:
            st.subheader("Two-Group Comparison Analysis (Volcano Plot + Bar Chart)")
            if detect_species_columns_from_diff is None:
                st.error("multivariate_stats 模块未正确加载")
            else:
                diff_df = st.session_state.get('diff_df')
                if diff_df is None:
                    st.warning("请先在侧边栏上传文件并点击 Start Analysis")
                else:
                    all_cols = diff_df.columns.tolist()
                    available_groups = []
                    for sp in ['LJ', 'BC', 'XY', 'EH', 'MB']:
                        if any(f'{sp}_Root' in c for c in all_cols):
                            available_groups.append(sp)
                    col_t2_1, col_t2_2 = st.columns([1, 1])
                    with col_t2_1:
                        group1 = st.selectbox("Group 1 (Numerator)", available_groups, index=0)
                    with col_t2_2:
                        group2 = st.selectbox("Group 2 (Denominator)", available_groups, index=1 if len(available_groups) > 1 else 0)
                    col_t2_3, col_t2_4, col_t2_5 = st.columns([1, 1, 1])
                    with col_t2_3:
                        p_thresh_2 = st.slider("P-value Threshold", 0.01, 0.10, 0.05, 0.01, key="p_thresh_2")
                    with col_t2_4:
                        vip_thresh_2 = st.slider("VIP Threshold", 0.5, 3.0, 1.0, 0.1, key="vip_thresh_2")
                    with col_t2_5:
                        fc_thresh_2 = st.slider("FC Threshold", 0.5, 3.0, 1.0, 0.1, key="fc_thresh_2")
                    search_metab = st.text_input("Search Metabolite (Highlight)", "")
                    chart_type = st.radio("Chart Type", ["Volcano Plot", "Bar Chart", "Both"], horizontal=True)
                    run_volcano_btn = st.button("Generate Comparison", type="primary")
                    if run_volcano_btn:
                        with st.spinner("Computing two-group comparison..."):
                            try:
                                cols_g1 = [c for c in all_cols if f'{group1}_Root' in c]
                                cols_g2 = [c for c in all_cols if f'{group2}_Root' in c]
                                if not cols_g1 or not cols_g2:
                                    st.error(f"未找到 {group1} 或 {group2} 的丰度列")
                                else:
                                    g1_mean = diff_df[cols_g1].mean(axis=1)
                                    g2_mean = diff_df[cols_g2].mean(axis=1)
                                    fc_vals = g1_mean / (g2_mean + 1e-9)
                                    diff_df['_fc_calc'] = fc_vals
                                    diff_df['_log2fc_calc'] = np.where(fc_vals >= 1, np.log2(fc_vals + 1e-9), -np.log2(1 / fc_vals + 1e-9))
                                    def classify_diff(row):
                                        vip = row.get('Vip_plsda', 0)
                                        pval = row.get('P_value', 1)
                                        fc = row.get('_fc_calc', 1)
                                        if vip > vip_thresh_2 and pval < p_thresh_2 and fc > fc_thresh_2:
                                            return 'Up'
                                        elif vip > vip_thresh_2 and pval < p_thresh_2 and fc < 1 / fc_thresh_2:
                                            return 'Down'
                                        return 'nosig'
                                    diff_df['_diff_class'] = diff_df.apply(classify_diff, axis=1)
                                    up_count = (diff_df['_diff_class'] == 'Up').sum()
                                    down_count = (diff_df['_diff_class'] == 'Down').sum()
                                    st.session_state['volcano_result'] = {'df': diff_df, 'group1': group1, 'group2': group2, 'p_thresh': p_thresh_2, 'vip_thresh': vip_thresh_2, 'fc_thresh': fc_thresh_2, 'up': up_count, 'down': down_count}
                                    st.success(f"Comparison: {group1} vs {group2} | Up: {up_count} | Down: {down_count}")
                            except Exception as e:
                                st.error(f"Comparison failed: {e}")
                    if 'volcano_result' in st.session_state:
                        res = st.session_state['volcano_result']
                        df_show = res['df'].copy()
                        if search_metab:
                            df_show['_highlight'] = df_show['Metabolite'].str.contains(search_metab, case=False, na=False)
                        else:
                            df_show['_highlight'] = False
                        if chart_type in ["Volcano Plot", "Both"]:
                            st.markdown(f"**Volcano Plot: {res['group1']} vs {res['group2']}**")
                            fig_vol = plot_volcano(df_show, p_thresh=res['p_thresh'], vip_thresh=res['vip_thresh'], fc_thresh=res['fc_thresh'])
                            st.plotly_chart(fig_vol, width="stretch")
                        if chart_type in ["Bar Chart", "Both"]:
                            st.markdown(f"**Differential Metabolite Counts: {res['group1']} vs {res['group2']}**")
                            fig_bar = plot_diff_bar(res['up'], res['down'], f"{res['group1']}_vs_{res['group2']}")
                            st.plotly_chart(fig_bar, width="stretch")
                        st.markdown("---")
                        st.subheader("Differential Metabolites Table")
                        sig_df = df_show[df_show['_diff_class'] != 'nosig'].copy().sort_values('Vip_plsda', ascending=False)
                        disp_cols = ['Metabolite', 'Vip_plsda', 'P_value', 'fdr', '_fc_calc', '_log2fc_calc', '_diff_class']
                        disp_cols = [c for c in disp_cols if c in sig_df.columns]
                        st.dataframe(sig_df[disp_cols], height=400)
                        if len(sig_df) > 0:
                            dl_df = sig_df[disp_cols].copy()
                            dl_bytes = df_to_excel_bytes(dl_df)
                            st.download_button(f"Download Differential Metabolites ({res['group1']}_vs_{res['group2']})", data=dl_bytes, file_name=f"differential_metabolites_{res['group1']}_vs_{res['group2']}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        with tab8:
            st.subheader("Multi-Group Comparison Analysis (ANOVA / Kruskal-Wallis)")
            if detect_species_columns_from_diff is None:
                st.error("multivariate_stats 模块未正确加载")
            else:
                diff_df = st.session_state.get('diff_df')
                if diff_df is None:
                    st.warning("请先在侧边栏上传文件并点击 Start Analysis")
                else:
                    all_cols = diff_df.columns.tolist()
                    available_groups = []
                    for sp in ['LJ', 'BC', 'XY', 'EH', 'MB']:
                        if any(f'{sp}_Root' in c for c in all_cols):
                            available_groups.append(sp)
                    col_t3_1, col_t3_2 = st.columns([1, 1])
                    with col_t3_1:
                        selected_groups_multi = st.multiselect("Select Groups for Comparison", available_groups, default=available_groups)
                    with col_t3_2:
                        test_method = st.selectbox("Statistical Test", ["Kruskal-Wallis H", "One-way ANOVA", "Both"])
                    col_t3_3, col_t3_4 = st.columns([1, 1])
                    with col_t3_3:
                        p_thresh_multi = st.slider("P-value Threshold", 0.01, 0.10, 0.05, 0.01, key="p_thresh_multi")
                    with col_t3_4:
                        top_n_heatmap = st.slider("Top N (Heatmap)", 20, 100, 50, 10, key="top_n_heatmap")
                    cluster_method = st.selectbox("Clustering Method", ["ward", "single", "complete", "average"])
                    run_multi_btn = st.button("Run Multi-Group Analysis", type="primary")
                    if run_multi_btn:
                        if len(selected_groups_multi) < 3:
                            st.warning("请至少选择3个组进行分析")
                        else:
                            with st.spinner("Running multi-group analysis..."):
                                try:
                                    groups_dict = {}
                                    for sp in selected_groups_multi:
                                        cols = [c for c in all_cols if f'{sp}_Root' in c]
                                        if cols:
                                            groups_dict[sp] = cols
                                    rows = []
                                    for sp, cols in groups_dict.items():
                                        for col in cols:
                                            row = {'Sample': col, 'Group': sp}
                                            for _, mrow in diff_df.iterrows():
                                                row[mrow['Metabolite']] = mrow[col]
                                            rows.append(row)
                                    abundance_df = pd.DataFrame(rows).set_index('Sample')
                                    metabolites = diff_df['Metabolite'].unique()
                                    test_results = []
                                    for metab in metabolites:
                                        if metab not in abundance_df.columns:
                                            continue
                                        groups_data = []
                                        valid_groups = []
                                        for sp in selected_groups_multi:
                                            if sp in groups_dict:
                                                cols = groups_dict[sp]
                                                vals = pd.to_numeric(abundance_df.loc[cols, metab], errors='coerce').dropna().values
                                                if len(vals) >= 2:
                                                    groups_data.append(vals)
                                                    valid_groups.append(sp)
                                        if len(groups_data) < 2:
                                            continue
                                        try:
                                            if test_method in ["Kruskal-Wallis H", "Both"]:
                                                h_stat, p_kw = stats.kruskal(*groups_data)
                                            else:
                                                h_stat, p_kw = np.nan, np.nan
                                            if test_method in ["One-way ANOVA", "Both"]:
                                                f_stat, p_anova = stats.f_oneway(*groups_data)
                                            else:
                                                f_stat, p_anova = np.nan, np.nan
                                            test_results.append({'Metabolite': metab, 'F_statistic': f_stat, 'P_ANOVA': p_anova, 'H_statistic': h_stat, 'P_Kruskal': p_kw, 'N_groups': len(valid_groups)})
                                        except Exception:
                                            continue
                                    results_df = pd.DataFrame(test_results)
                                    if test_method in ["Kruskal-Wallis H", "Both"] and 'P_Kruskal' in results_df.columns:
                                        results_df = results_df.sort_values('P_Kruskal')
                                        n = len(results_df)
                                        results_df['P_Kruskal_BH'] = results_df['P_Kruskal'] * n / (results_df.index + 1)
                                    if test_method in ["One-way ANOVA", "Both"] and 'P_ANOVA' in results_df.columns:
                                        results_df = results_df.sort_values('P_ANOVA')
                                        n = len(results_df)
                                        results_df['P_ANOVA_BH'] = results_df['P_ANOVA'] * n / (results_df.index + 1)
                                    vip_map = diff_df.groupby('Metabolite')['Vip_plsda'].mean().to_dict()
                                    results_df['VIP_mean'] = results_df['Metabolite'].map(vip_map)
                                    st.session_state['multi_group_results'] = {'results': results_df, 'groups_dict': groups_dict, 'p_thresh': p_thresh_multi, 'top_n': top_n_heatmap, 'cluster_method': cluster_method}
                                    sig_count = (results_df['P_Kruskal'] < p_thresh_multi).sum() if 'P_Kruskal' in results_df.columns else 0
                                    st.success(f"Multi-group analysis completed! Significant metabolites: {sig_count} (P < {p_thresh_multi})")
                                except Exception as e:
                                    st.error(f"Analysis failed: {e}")
                    if 'multi_group_results' in st.session_state:
                        res = st.session_state['multi_group_results']
                        results_df = res['results']
                        if 'P_Kruskal' in results_df.columns:
                            sig_df = results_df[results_df['P_Kruskal'] < res['p_thresh']].copy()
                        elif 'P_ANOVA' in results_df.columns:
                            sig_df = results_df[results_df['P_ANOVA'] < res['p_thresh']].copy()
                        else:
                            sig_df = results_df.copy()
                        sig_df = sig_df.sort_values('P_Kruskal' if 'P_Kruskal' in sig_df.columns else 'P_ANOVA')
                        st.markdown(f"**Significant Metabolites: {len(sig_df)}**")
                        disp_cols = ['Metabolite', 'F_statistic', 'P_ANOVA', 'H_statistic', 'P_Kruskal', 'VIP_mean']
                        disp_cols = [c for c in disp_cols if c in sig_df.columns]
                        st.dataframe(sig_df[disp_cols].head(50), height=400)
                        st.markdown("---")
                        st.subheader("Clustered Heatmap")
                        top_metabs = sig_df.head(res['top_n'])['Metabolite'].tolist() if len(sig_df) > 0 else results_df.head(res['top_n'])['Metabolite'].tolist()
                        if top_metabs:
                            try:
                                all_samples = []
                                for sp, cols in res['groups_dict'].items():
                                    all_samples.extend(cols)
                                heatmap_data = diff_df[diff_df['Metabolite'].isin(top_metabs)][['Metabolite'] + all_samples].copy().set_index('Metabolite')
                                z_data = heatmap_data.apply(lambda x: (x - x.mean()) / (x.std() + 1e-9), axis=1).fillna(0)
                                z_data.index = [n[:25] + '...' if len(n) > 25 else n for n in z_data.index]
                                fig_heat = px.imshow(z_data.values, x=z_data.columns, y=z_data.index, color_continuous_scale='RdBu_r', title=f'Clustered Heatmap (Top {len(top_metabs)} by P-value, Z-score)', labels=dict(x='Sample', y='Metabolite', color='Z-score'), height=max(400, len(z_data) * 8), width=900)
                                fig_heat.update_layout(xaxis={'tickangle': 45}, font=dict(size=9))
                                st.plotly_chart(fig_heat, width="stretch")
                            except Exception as e:
                                st.warning(f"Heatmap generation failed: {e}")
                        st.markdown("---")
                        st.subheader("Boxplot of Top Differential Metabolites")
                        if len(sig_df) > 0:
                            top5_metabs = sig_df.head(5)['Metabolite'].tolist()
                            metab_to_show = st.selectbox("Select Metabolite", top5_metabs)
                            try:
                                rows_box = []
                                for sp, cols in res['groups_dict'].items():
                                    for col in cols:
                                        val = diff_df[diff_df['Metabolite'] == metab_to_show][col].values
                                        if len(val) > 0:
                                            rows_box.append({'Group': sp, 'Sample': col, 'Abundance': val[0]})
                                box_df = pd.DataFrame(rows_box)
                                box_df['Abundance'] = pd.to_numeric(box_df['Abundance'], errors='coerce')
                                box_df = box_df.dropna()
                                if len(box_df) > 0:
                                    fig_box = px.box(box_df, x='Group', y='Abundance', color='Group', title=f'Abundance Distribution: {metab_to_show[:40]}', color_discrete_map={g: GROUP_COLORS.get(g, '#3498db') for g in box_df['Group'].unique()})
                                    fig_box.update_layout(height=400)
                                    st.plotly_chart(fig_box, width="stretch")
                            except Exception as e:
                                st.warning(f"Boxplot generation failed: {e}")
                        if len(results_df) > 0:
                            dl_bytes = df_to_excel_bytes(results_df.copy())
                            st.download_button("Download Multi-Group Analysis Results", data=dl_bytes, file_name="multi_group_comparison_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        with tab9:
            st.subheader("🧬 Network Pharmacology — SMILES / Target Prediction / Enrichment")
            if query_smiles_by_name is None:
                st.error("network_pharma 模块未正确加载，请检查 network_pharma.py 是否与 app.py 在同一目录。")
            else:
                # ------- 子模块选择 -------
                np_mode = st.radio(
                    "选择子模块",
                    ["① SMILES 查询", "② SwissTargetPrediction 靶点预测", "③ gseapy 富集分析"],
                    horizontal=True,
                )

                # ====== ① SMILES 查询 ======
                if "①" in np_mode:
                    st.markdown("#### ① 通过化合物名称查询 SMILES（PubChem）")
                    st.caption("输入化合物名称或 CAS 号，自动查询 SMILES 结构式。结果缓存 1 小时。")
                    sm_name = st.text_input("化合物名称或 CAS 号", placeholder="例如：Quercetin 或 117-39-5")
                    if st.button("查询 SMILES", type="primary"):
                        if not sm_name.strip():
                            st.warning("请输入化合物名称")
                        else:
                            with st.spinner("查询中..."):
                                sm = query_smiles_by_name(sm_name.strip())
                            if sm != "NOT_FOUND":
                                st.success(f"**SMILES:** `{sm}`")
                                st.code(sm, language="text")
                            else:
                                st.warning("未找到该化合物的 SMILES，请尝试其他名称或检查拼写")

                    st.markdown("---")
                    st.markdown("#### 批量 SMILES 查询（基于当前候选化合物）")
                    if st.button("为所有药理候选化合物查询 SMILES", type="secondary"):
                        if st.session_state.get('pharma_df') is None:
                            st.warning("请先在侧边栏运行分析")
                        else:
                            df_in = st.session_state['pharma_df'].copy()
                            df_in = df_in.rename(columns={'compound_name': 'compound_name'})
                            progress_bar = st.progress(0)
                            def pg_cb(cur, tot):
                                progress_bar.progress(int(cur / tot * 100))
                            df_out = batch_query_smiles(df_in, name_col='compound_name', progress_callback=pg_cb)
                            progress_bar.empty()
                            sm_found = (df_out['SMILES'] != 'NOT_FOUND').sum()
                            st.success(f"查询完成！找到 SMILES: {sm_found}/{len(df_out)}")
                            # 存入 session_state 供下游使用
                            st.session_state['smiles_df'] = df_out
                            st.dataframe(df_out[['Metabolite', 'compound_name', 'SMILES']].head(30), height=300)
                            b = io.BytesIO()
                            with pd.ExcelWriter(b, engine='openpyxl') as w:
                                df_out[['Metabolite', 'compound_name', 'SMILES']].to_excel(w, index=False, sheet_name='SMILES')
                            b.seek(0)
                            st.download_button("Download SMILES Results", data=b.getvalue(),
                                file_name="smiles_query_results.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

                # ====== ② SwissTargetPrediction ======
                elif "②" in np_mode:
                    st.markdown("#### ② SwissTargetPrediction 靶点预测")
                    st.caption("输入 SMILES，自动预测人类靶点基因（Probability > 0.01）。每个 SMILES 约需 1-2 秒。")
                    st_smiles = st.text_area("输入 SMILES（支持批量，每行一个）",
                        placeholder="CC1=CC=C(C=C1)C2=CC(=NN2C=3C=C(C=C3)C4=CC=CC=C4)C5=CC=CC=C5\nCC(O)=O")
                    if st.button("预测靶点", type="primary"):
                        if not st_smiles.strip():
                            st.warning("请输入至少一个 SMILES")
                        else:
                            smiles_list = [s.strip() for s in st_smiles.strip().split('\n') if s.strip()]
                            results = []
                            bar = st.progress(0)
                            for i, sm in enumerate(smiles_list):
                                genes, cnt, err = query_swiss_target_prediction(sm)
                                results.append({'SMILES': sm, 'Predicted_Targets': genes,
                                               'Target_Count': cnt, 'Error': err})
                                bar.progress(int((i+1)/len(smiles_list)*100))
                            bar.empty()
                            res_df = pd.DataFrame(results)
                            st.session_state['target_df'] = res_df
                            st.success(f"完成！{len(results)} 个化合物中有靶点记录: {(res_df['Target_Count']>0).sum()}")
                            st.dataframe(res_df, height=300)
                            b = io.BytesIO()
                            with pd.ExcelWriter(b, engine='openpyxl') as w:
                                res_df.to_excel(w, index=False, sheet_name='Targets')
                            b.seek(0)
                            st.download_button("Download Target Prediction", data=b.getvalue(),
                                file_name="swiss_target_prediction.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

                    st.markdown("---")
                    st.markdown("#### 批量预测（基于 SMILES 查询结果）")
                    if st.button("对 SMILES 查询结果批量预测靶点", type="secondary"):
                        if st.session_state.get('smiles_df') is None:
                            st.warning("请先运行 ① SMILES 查询生成数据")
                        else:
                            df_in = st.session_state['smiles_df'].copy()
                            # 过滤有效 SMILES
                            df_in = df_in[df_in['SMILES'].notna() & (df_in['SMILES'] != 'NOT_FOUND') & (df_in['SMILES'] != 'PENDING')]
                            if len(df_in) == 0:
                                st.warning("没有有效的 SMILES 可供预测")
                            else:
                                progress_bar2 = st.progress(0)
                                def pg_cb2(cur, tot):
                                    progress_bar2.progress(int(cur/tot*100))
                                df_out, errs = batch_swiss_target_prediction(df_in, smiles_col='SMILES',
                                                                           progress_callback=pg_cb2)
                                progress_bar2.empty()
                                n_hit = (df_out['Target_Count'] > 0).sum()
                                st.success(f"完成！{n_hit}/{len(df_out)} 个化合物命中靶点")
                                if errs:
                                    st.warning(f"出错 {len(errs)} 个: {errs[:3]}")
                                st.session_state['target_df'] = df_out
                                disp = df_out[df_out['Target_Count'] > 0][['Metabolite', 'SMILES', 'Predicted_Targets', 'Target_Count']]
                                st.dataframe(disp, height=300)
                                b = io.BytesIO()
                                with pd.ExcelWriter(b, engine='openpyxl') as w:
                                    df_out.to_excel(w, index=False, sheet_name='Targets')
                                b.seek(0)
                                st.download_button("Download Batch Targets", data=b.getvalue(),
                                    file_name="swiss_target_batch.xlsx",
                                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

                # ====== ③ gseapy 富集分析 ======
                elif "③" in np_mode:
                    st.markdown("#### ③ gseapy 富集分析（KEGG / GO）")
                    st.caption("输入基因列表，自动在 Enrichr 数据库做 KEGG + GO 富集分析（需要网络连接）。")
                    gene_input = st.text_area("输入基因列表（每行一个，或用分号/逗号分隔）",
                        placeholder="EGFR\nAKT1\nTP53\nVEGFA\nIL6\nTNF")
                    col_db1, col_db2 = st.columns(2)
                    with col_db1:
                        organisms = st.selectbox("物种", ["Human", "Mouse", "Rat"], index=0)
                    with col_db2:
                        top_n_enr = st.number_input("每库显示 Top N", 5, 50, 20)
                    if st.button("运行富集分析", type="primary"):
                        if not gene_input.strip():
                            st.warning("请输入基因列表")
                        else:
                            # 解析基因列表
                            raw = gene_input.replace(';', '\n').replace(',', '\n')
                            genes = [g.strip() for g in raw.split('\n') if g.strip()]
                            if not genes:
                                st.warning("未能解析有效基因")
                            else:
                                org_map = {"Human": "human", "Mouse": "mouse", "Rat": "rat"}
                                with st.spinner(f"正在分析 {len(genes)} 个基因..."):
                                    try:
                                        enr_dict = run_gseapy_enrichment(
                                            gene_list=genes,
                                            organism=org_map[organisms],
                                        )
                                        merged = merge_enrichment_results(enr_dict)
                                        st.session_state['enrichment_df'] = merged
                                        if len(merged) == 0:
                                            st.warning("未找到显著富集结果")
                                        else:
                                            st.success(f"富集分析完成，共 {len(merged)} 条结果")
                                            # 汇总统计
                                            for db, grp in merged.groupby('Database'):
                                                st.markdown(f"**{db}:** {len(grp)} 条结果")
                                            # Dotplot
                                            if plot_enrichment_dotplot is not None:
                                                fig = plot_enrichment_dotplot(merged.head(top_n_enr * 2))
                                                if fig:
                                                    st.plotly_chart(fig, width="stretch")
                                            st.dataframe(merged[['Database', 'Term', 'P_value', 'Adjusted_P_value', 'Genes', 'Overlap']].head(50),
                                                         height=400)
                                            b = io.BytesIO()
                                            with pd.ExcelWriter(b, engine='openpyxl') as w:
                                                merged.to_excel(w, index=False, sheet_name='Enrichment')
                                            b.seek(0)
                                            st.download_button("Download Enrichment Results",
                                                data=b.getvalue(),
                                                file_name="gseapy_enrichment.xlsx",
                                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                                    except Exception as e:
                                        st.error(f"富集分析失败: {e}")

                    # ---- 从靶点预测结果直接导入基因列表 ----
                    st.markdown("---")
                    st.markdown("#### 一键分析：使用 SwissTargetPrediction 靶点基因")
                    if st.button("用当前靶点预测结果做富集分析", type="secondary"):
                        if st.session_state.get('target_df') is None:
                            st.warning("请先运行 ② SwissTargetPrediction 生成靶点数据")
                        else:
                            # 合并所有靶点基因
                            all_genes = set()
                            for _, row in st.session_state['target_df'].iterrows():
                                tg = str(row.get('Predicted_Targets', ''))
                                if tg and tg not in ('NOT_FOUND', 'ERROR', ''):
                                    for g in tg.split(';'):
                                        g = g.strip()
                                        if g:
                                            all_genes.add(g)
                            if not all_genes:
                                st.warning("没有找到有效靶点基因")
                            else:
                                gene_list = sorted(all_genes)
                                st.info(f"汇总了 {len(gene_list)} 个靶点基因，开始富集分析...")
                                org_map = {"Human": "human", "Mouse": "mouse", "Rat": "rat"}
                                with st.spinner(f"正在分析 {len(gene_list)} 个基因..."):
                                    try:
                                        enr_dict = run_gseapy_enrichment(gene_list=gene_list, organism="human")
                                        merged = merge_enrichment_results(enr_dict)
                                        st.session_state['enrichment_df'] = merged
                                        if len(merged) == 0:
                                            st.warning("未找到显著富集结果")
                                        else:
                                            st.success(f"富集分析完成，共 {len(merged)} 条结果（4个数据库）")
                                            if plot_enrichment_dotplot is not None:
                                                fig = plot_enrichment_dotplot(merged.head(top_n_enr * 2))
                                                if fig:
                                                    st.plotly_chart(fig, width="stretch")
                                            st.dataframe(merged[['Database', 'Term', 'P_value', 'Adjusted_P_value', 'Genes', 'Overlap']].head(50),
                                                         height=400)
                                            b = io.BytesIO()
                                            with pd.ExcelWriter(b, engine='openpyxl') as w:
                                                merged.to_excel(w, index=False, sheet_name='Enrichment')
                                            b.seek(0)
                                            st.download_button("Download Enrichment Results",
                                                data=b.getvalue(),
                                                file_name="gseapy_enrichment.xlsx",
                                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                                    except Exception as e:
                                        st.error(f"富集分析失败: {e}")

    else:
        st.info("Please upload files in the sidebar and click **Start Analysis** to begin.")
        st.markdown("---")
        st.subheader("Workflow")
        st.markdown("""
        **Step 1: Upload Files** - Upload all `*_diff.exp.xls` comparison files (multiple files supported) + `anno.xls` annotation file
        **Step 2: Configure Parameters** - Adjust VIP and P-value thresholds in the sidebar
        **Step 3: Analysis** - Click **Start Analysis** to run the pipeline
        **Step 4: Results** - 📊 Data Overview · 📈 Visualization · ⭐ Star Molecules · 🧬 KEGG Enrichment
        """)


if __name__ == '__main__':
    main()
