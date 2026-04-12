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
        load_tcmsp, load_drugbank, load_ttd,
        match_pharma_db, compute_pharma_evidence_score
    )
except ImportError:
    load_tcmsp = load_drugbank = load_ttd = match_pharma_db = compute_pharma_evidence_score = None
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

# 中文字体支持
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
    # ===== 侧边栏 =====
    st.sidebar.title("Config")
    st.sidebar.markdown("---")

    # 文件上传
    st.sidebar.subheader("Upload Files")
    uploaded_diff = st.sidebar.file_uploader(
        "上传 diff.exp.xls 文件（支持多选）",
        type=['xls', 'xlsx', 'tsv', 'txt'],
        accept_multiple_files=True,
        help="同时选中所有比较组的文件，如 BC_vs_LJ.diff.exp.xls、LJ_vs_XY.diff.exp.xls 等"
    )
    uploaded_anno = st.sidebar.file_uploader(
        "上传 anno.xls 注释文件",
        type=['xls', 'xlsx', 'tsv', 'txt'],
        help="包含化合物注释信息的anno.xls文件"
    )

    st.sidebar.markdown("---")
    st.sidebar.subheader("Parameters")

    vip_thresh = st.sidebar.slider(
        "VIP Threshold (Vip_plsda > ?)", 0.5, 3.0, 1.0, 0.1,
        help="只保留 VIP > 此值的差异代谢物"
    )
    p_thresh = st.sidebar.slider(
        "P-value Threshold (P_value < ?)", 0.01, 0.1, 0.05, 0.005,
        help="只保留 P_value < 此值的差异代谢物"
    )

    run_analysis = st.sidebar.button("Start Analysis", type="primary", use_container_width=True)

    st.sidebar.markdown("---")
    st.sidebar.markdown(
        "**Input format**: `.xls`/`.xlsx`/`.tsv` 文件（自动识别）"
    )

    # ===== 主区域：标题 =====
    st.title("差异代谢物药用筛选 Web")
    st.markdown("上传差异分析文件 + 注释文件，自动完成筛选、可视化、明星分子精选、KEGG 富集分析。")

    # ===== 状态初始化 =====
    if 'analysis_done' not in st.session_state:
        st.session_state['analysis_done'] = False
        st.session_state['diff_df'] = None
        st.session_state['pharma_df'] = None
        st.session_state['abundance_df'] = None
        st.session_state['annotated_df'] = None
        st.session_state['star_df'] = None
        st.session_state['uploaded_diff'] = None
        st.session_state['uploaded_anno'] = None

    # ===== 运行分析 =====
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

        # 保存到 session state
        st.session_state['diff_df'] = merged
        st.session_state['diff_filtered'] = diff_filtered
        st.session_state['annotated_df'] = annotated
        st.session_state['pharma_df'] = pharma
        st.session_state['abundance_df'] = abundance
        st.session_state['star_df'] = star
        st.session_state['analysis_done'] = True
        st.session_state['uploaded_diff'] = uploaded_diff
        st.session_state['uploaded_anno'] = uploaded_anno
        st.success("Analysis complete!")

    # ===== 主区域：分析结果展示 =====
    if st.session_state['analysis_done']:
        # ===== 标签页 =====
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "📊 Data Overview & Downloads",
            "📈 Visualization",
            "⭐ Star Molecules",
            "🧬 KEGG Enrichment",
            "💊 Pharma DB",
        ])

        # ---------- Tab 1: 数据概览 + 下载 ----------
        with tab1:
            st.subheader("Data Overview")

            col_stat1, col_stat2, col_stat3, col_stat4 = st.columns(4)
            col_stat1.metric("Total Metabolites", len(st.session_state['diff_df']))
            col_stat2.metric("Differential (VIP>P)", len(st.session_state['diff_filtered']))
            named = st.session_state['annotated_df']['compound_name'].notna() & \
                    (st.session_state['annotated_df']['compound_name'] != '-')
            col_stat3.metric("Named Compounds", int(named.sum()))
            col_stat4.metric("Pharmacological", len(st.session_state['pharma_df']))

            st.markdown("---")

            # 三个下载按钮
            st.subheader("Download Results")

            col_dl1, col_dl2, col_dl3 = st.columns(3)

            # 下载1: 总差异代谢物
            diff_out_cols = ['Metabolite', 'Vip_plsda', 'Vip_oplsda', 'P_value', 'fdr', 'FC',
                             'compound_name', 'formula', 'super_class', 'kegg_second_category']
            diff_out_cols = [c for c in diff_out_cols if c in st.session_state['annotated_df'].columns]
            diff_out = st.session_state['annotated_df'][diff_out_cols].drop_duplicates(subset=['Metabolite'])
            bytes1 = df_to_excel_bytes(diff_out)

            col_dl1.download_button(
                label="Download: Total Differential Metabolites",
                data=bytes1,
                file_name="total_differential_metabolites.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
            )
            col_dl1.caption(f"{len(diff_out)} rows")

            # 下载2: 候选药用化合物
            out2_cols = ['Metabolite', 'compound_name', 'formula', 'hmdb_id',
                         'super_class', 'class', 'sub_class',
                         'compound_first_category', 'compound_second_category',
                         'kegg_first_category', 'kegg_second_category',
                         'pathway_id', 'Vip_plsda', 'P_value', 'FC', '_pharma_reason']
            out2_cols = [c for c in out2_cols if c in st.session_state['pharma_df'].columns]
            pharma_out = st.session_state['pharma_df'][out2_cols].copy()
            pharma_out = pharma_out.rename(columns={'_pharma_reason': '筛选依据'})
            bytes2 = df_to_excel_bytes(pharma_out)

            col_dl2.download_button(
                label="Download: Pharmacological Candidates",
                data=bytes2,
                file_name="pharmacological_candidates.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
            )
            col_dl2.caption(f"{len(pharma_out)} rows")

            # 下载3: 五物种丰度
            bytes3 = df_to_excel_bytes(st.session_state['abundance_df'])
            col_dl3.download_button(
                label="Download: Five-Species Abundance",
                data=bytes3,
                file_name="five_species_abundance.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
            )
            col_dl3.caption(f"{len(st.session_state['abundance_df'])} rows")

            st.markdown("---")

            # ===== Word 报告下载 =====
            if generate_word_report is not None:
                st.subheader("Download Analysis Report")
                if st.button("Generate Word Report", type="secondary", use_container_width=True):
                    with st.spinner("Generating Word report..."):
                        try:
                            fig_hm_bytes = None
                            fig_bx_bytes = None
                            fig_vip_bytes = None
                            fig_bubble_bytes = None
                            fig_bar_bytes = None
                            # 尝试复用已生成的图表字节
                            try:
                                from io import BytesIO
                                buf = BytesIO()
                                pharma_df = st.session_state['pharma_df']
                                abundance_df = st.session_state['abundance_df']
                                annotated_df = st.session_state['annotated_df']
                                star_df = st.session_state.get('star_df', pharma_df)

                                params = {
                                    'vip_thresh': st.session_state.get('vip_thresh', 1.0),
                                    'p_thresh': st.session_state.get('p_thresh', 0.05),
                                    'n_diff': len(st.session_state['diff_filtered']),
                                    'n_named': int((annotated_df['compound_name'].notna() & (annotated_df['compound_name'] != '-')).sum()),
                                    'n_pharma': len(pharma_df),
                                    'n_star': len(star_df),
                                }
                                doc_bytes = generate_word_report(
                                    annotated_df=annotated_df,
                                    pharma_df=pharma_df,
                                    star_df=star_df,
                                    abundance_df=abundance_df,
                                    params=params,
                                    fig_heatmap_bytes=fig_hm_bytes,
                                    fig_boxplot_bytes=fig_bx_bytes,
                                    fig_vip_bytes=fig_vip_bytes,
                                    fig_bubble_bytes=fig_bubble_bytes,
                                    fig_enr_bar_bytes=fig_bar_bytes,
                                    enr_df=st.session_state.get('enr_df'),
                                )
                            except Exception as e:
                                doc_bytes = None
                                st.warning(f"Report generation partially failed: {e}")

                        except Exception as e:
                            doc_bytes = None
                            st.error(f"Report generation failed: {e}")

                    if doc_bytes:
                        st.success("Report generated successfully!")
                        st.download_button(
                            "Download Word Report (.docx)",
                            data=doc_bytes,
                            file_name=f"metabolite_analysis_report_{datetime.now().strftime('%Y%m%d_%H%M')}.docx",
                            mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                            use_container_width=True,
                        )
            else:
                st.caption("Note: python-docx not installed. Run: pip install python-docx")

            st.markdown("---")
            st.subheader("Preview: Pharmacological Candidates")
            display_cols = ['Metabolite', 'compound_name', 'super_class', 'Vip_plsda', 'P_value', 'FC', '_pharma_reason']
            display_cols = [c for c in display_cols if c in st.session_state['pharma_df'].columns]
            st.dataframe(
                st.session_state['pharma_df'][display_cols].head(30),
                use_container_width=True,
                height=400,
            )

        # ---------- Tab 2: 可视化 ----------
        with tab2:
            st.subheader("Visualization")

            viz_col1, viz_col2 = st.columns([1, 1])

            # 热图
            with viz_col1:
                st.markdown("**Clustered Heatmap (Top N by VIP)**")
                hm_top_n = st.slider("Select Top N compounds for heatmap:", 10, 100, 30, 5, key="hm_top")
                fig_hm = plot_heatmap(st.session_state['pharma_df'], st.session_state['abundance_df'], top_n=hm_top_n)
                if fig_hm:
                    st.pyplot(fig_hm)
                    buf_fig = io.BytesIO()
                    fig_hm.savefig(buf_fig, format='png', dpi=150, bbox_inches='tight')
                    buf_fig.seek(0)
                    st.download_button(
                        "Download Heatmap PNG",
                        data=buf_fig,
                        file_name="heatmap.png",
                        mime="image/png",
                        use_container_width=True,
                    )
                else:
                    st.warning("Insufficient data for heatmap.")

            # 箱线图
            with viz_col2:
                st.markdown("**Boxplot by Species (Top N by VIP)**")
                bx_top_n = st.slider("Select Top N compounds for boxplot:", 5, 50, 15, 5, key="bx_top")
                fig_bx = plot_boxplot(st.session_state['abundance_df'], st.session_state['pharma_df'], top_n=bx_top_n)
                if fig_bx:
                    st.plotly_chart(fig_bx, use_container_width=True)
                else:
                    st.warning("Insufficient data for boxplot.")

            st.markdown("---")

            # VIP 分布图
            st.subheader("VIP Score Distribution")
            fig_vip = px.histogram(
                st.session_state['annotated_df'],
                x='Vip_plsda',
                nbins=50,
                title='VIP Score Distribution (Differential Metabolites)',
                labels={'Vip_plsda': 'VIP (PLS-DA)', 'count': 'Count'},
                color_discrete_sequence=['#4C78A8'],
            )
            fig_vip.add_vline(x=vip_thresh, line_dash='dash', line_color='red',
                              annotation_text=f'VIP threshold={vip_thresh}')
            fig_vip.update_layout(height=400)
            st.plotly_chart(fig_vip, use_container_width=True)

        # ---------- Tab 3: 明星分子精选 ----------
        with tab3:
            st.subheader("Star Molecule Selection")

            # 优先使用含药理证据的 v2 评分，否则用原始评分
            if 'star_df_v2' in st.session_state:
                star_df = st.session_state['star_df_v2'].copy()
                score_col = '_star_score_v2'
                st.info("已更新为含药理证据的综合评分 (v2)")
            else:
                star_df = st.session_state['star_df'].copy()
                score_col = '_star_score'

            all_super_classes = sorted(star_df['super_class'].dropna().unique().tolist())

            # 筛选控件
            col_f1, col_f2, col_f3 = st.columns([1, 1, 1])

            with col_f1:
                selected_classes = st.multiselect(
                    "Pharmacological Category",
                    options=all_super_classes,
                    default=['Alkaloids and derivatives', 'Phenylpropanoids and polyketides'],
                    format_func=lambda x: x if x else 'All',
                    help="Select categories to include"
                )

            with col_f2:
                min_vip = st.slider(
                    "Minimum VIP", 0.5, 5.0, 1.0, 0.1,
                    help="Only show compounds with VIP above this value"
                )

            with col_f3:
                min_diff_ratio = st.slider(
                    "Minimum Species Diff Ratio (max/min)", 1.0, 100.0, 1.0, 0.5,
                    help="Filter by maximum/minimum abundance ratio"
                )

            # 应用筛选
            filtered_star = star_df[star_df['Vip_plsda'] >= min_vip].copy()

            if selected_classes:
                filtered_star = filtered_star[filtered_star['super_class'].isin(selected_classes)]

            if '_diff_ratio' in filtered_star.columns:
                filtered_star = filtered_star[filtered_star['_diff_ratio'] >= min_diff_ratio]

            st.markdown(f"Filtered: **{len(filtered_star)}** / {len(star_df)} candidates")

            # Top 10 明星分子
            st.markdown("---")
            st.subheader("Top 10 Star Molecules (by Composite Score)")

            top10 = filtered_star.head(10)
            cols_show = ['Metabolite', 'compound_name', 'super_class', 'Vip_plsda',
                         '_diff_ratio', score_col, '_pharma_reason']
            cols_show = [c for c in cols_show if c in top10.columns]

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
                        st.markdown(
                            f"**#{rank} {metab_name}**  "
                            f"| VIP={vip:.2f} | Diff={diff_r:.1f}x | Score={score:.3f}  "
                            f"| *{sc}*  "
                            f"| {reason}"
                        )
                    with c2:
                        # 文献链接
                        search_name = str(compound) if compound and compound != '-' else str(metab_name)
                        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={search_name.replace(' ', '+')}"
                        scholar_url = f"https://scholar.google.com/scholar?q={search_name.replace(' ', '+')}"
                        st.markdown(
                            f"[PubMed]({pubmed_url}) · "
                            f"[Google Scholar]({scholar_url})",
                            unsafe_allow_html=True
                        )
                    st.divider()

            st.markdown("---")

            # 详细表格
            st.subheader("Full Filtered List (with Scores)")
            table_cols = ['Metabolite', 'compound_name', 'super_class', 'Vip_plsda',
                          '_diff_ratio', score_col, '_vip_norm', '_diff_norm', '_pharma_score',
                          'LJ_mean_abundance', 'BC_mean_abundance', 'XY_mean_abundance',
                          'EH_mean_abundance', 'MB_mean_abundance']
            table_cols = [c for c in table_cols if c in filtered_star.columns]
            st.dataframe(
                filtered_star[table_cols].style.format({
                    'Vip_plsda': '{:.4f}',
                    '_diff_ratio': '{:.2f}',
                    score_col: '{:.3f}',
                    '_vip_norm': '{:.3f}',
                    '_diff_norm': '{:.3f}',
                    '_pharma_score': '{:.2f}',
                    'LJ_mean_abundance': '{:.4f}',
                    'BC_mean_abundance': '{:.4f}',
                    'XY_mean_abundance': '{:.4f}',
                    'EH_mean_abundance': '{:.4f}',
                    'MB_mean_abundance': '{:.4f}',
                }, na_rep='-'),
                use_container_width=True,
                height=500,
            )

            # 下载精选结果
            bytes_star = df_to_excel_bytes(filtered_star[table_cols])
            st.download_button(
                "Download Star Molecules Table",
                data=bytes_star,
                file_name="star_molecules.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
            )

        # ---------- Tab 4: KEGG 通路富集 ----------
        with tab4:
            st.subheader("KEGG Pathway Enrichment Analysis")
            st.markdown(
                "基于注释数据计算 KEGG 通路过representation分析（超几何检验），"
                "无需联网。并提供 KEGG Mapper 在线可视化链接。"
            )
            st.info(
                "提示：KEGG Mapper 和 MetaboAnalyst 在线分析可获得通路图高亮。\n"
                "下方表格中点击 KEGG ID 可跳转到 KEGG Mapper 进行通路可视化。"
            )

            if st.button("Run Pathway Enrichment", type="primary", use_container_width=True):
                pharma_df = st.session_state['pharma_df']
                annotated_df = st.session_state['annotated_df']

                # Step 1: 构建代谢物->通路映射（从背景注释表）
                metab_to_pathways = {}
                for _, row in annotated_df.iterrows():
                    metab = row['Metabolite']
                    pid = str(row.get('pathway_id', '-'))
                    if pid not in ('-', 'nan', '') and pd.notna(row.get('pathway_id')):
                        pathways = [p.strip() for p in pid.split(';') if p.strip().startswith('map')]
                        if pathways:
                            metab_to_pathways[metab] = pathways

                # 背景：所有注释到KEGG通路的差异代谢物
                bg_mets = set(metab_to_pathways.keys())
                N_bg = len(bg_mets)

                # 统计每个通路在背景中的代谢物数
                pathway_bg_count = {}
                for pathways in metab_to_pathways.values():
                    for p in pathways:
                        pathway_bg_count[p] = pathway_bg_count.get(p, 0) + 1

                # Step 2: 候选化合物集合
                pharma_mets = set(pharma_df['Metabolite'].unique())
                n_cand = len(pharma_mets & set(metab_to_pathways.keys()))

                st.info(
                    f"背景代谢物（注释到KEGG通路）: {N_bg} | "
                    f"候选化合物（注释到KEGG通路）: {n_cand}"
                )

                # Step 3: 超几何检验（代谢物视角）
                # N = 背景代谢物总数（注释到KEGG），K = 该通路在背景中的代谢物数
                # n = 候选代谢物数，k = 候选中该通路命中数
                enr_results = []
                for pathway_id, K in pathway_bg_count.items():
                    hit_mets = [
                        m for m in pharma_mets
                        if m in metab_to_pathways and pathway_id in metab_to_pathways[m]
                    ]
                    k = len(hit_mets)
                    if k > 0:
                        pval = stats.hypergeom.sf(k - 1, N_bg, K, n_cand)
                        hit_names = '; '.join(sorted(hit_mets)[:8])
                        if len(hit_mets) > 8:
                            hit_names += f' ... (+{len(hit_mets) - 8} more)'

                        enr_results.append({
                            'KEGG_ID': pathway_id,
                            'KEGG_URL': f'https://www.genome.jp/kegg/mapper/?org_name=hsa&mode=compound&pathway={pathway_id}',
                            'Pathway_Name': _PATHWAY_NAME_MAP.get(pathway_id, pathway_id),
                            'Candidate_Hits': k,
                            'Total_in_Pathway': K,
                            'Candidate_Total': n_cand,
                            'Rich_Factor': round(k / K, 3) if K > 0 else 0,
                            'P_value': pval,
                            'Hit_Metabolites': hit_names,
                    })

                enr_df = pd.DataFrame(enr_results)
                if len(enr_df) == 0:
                    st.warning(
                        "候选化合物中未找到带 KEGG pathway_id 的注释，或无任何通路被命中。"
                        "请检查注释文件是否包含有效的 pathway_id。"
                    )
                elif len(enr_df) > 0:
                    # BH 校正 (Benjamini-Hochberg)
                    enr_df = enr_df.sort_values('P_value').reset_index(drop=True)
                    n_tests = len(enr_df)
                    # 按 P-value 排序后，i 即为 rank（1-indexed）
                    adj_pvals = []
                    for i, (_, row) in enumerate(enr_df.iterrows()):
                        rank = i + 1
                        adj_p = min(row['P_value'] * n_tests / rank, 1.0)
                        adj_pvals.append(adj_p)
                    enr_df['P_value_BH'] = adj_pvals
                    enr_df['-log10P'] = -np.log10(enr_df['P_value'].replace(0, 1e-10))

                    # 过滤 P < 0.05
                    sig_df = enr_df[enr_df['P_value'] < 0.05].copy()
                    st.success(f"显著富集通路: {len(sig_df)} 条 (P < 0.05)")

                    # ---- 气泡图 ----
                    st.markdown("**Pathway Enrichment Bubble Chart (P < 0.05)**")
                    plot_df = sig_df.head(30).copy()
                    plot_df['label'] = plot_df['KEGG_ID'] + ': ' + plot_df['Pathway_Name']
                    plot_df['label'] = plot_df['label'].apply(
                        lambda x: x[:55] + '...' if len(str(x)) > 55 else x
                    )
                    fig_bubble = px.scatter(
                        plot_df,
                        x='Rich_Factor',
                        y='label',
                        size='Candidate_Hits',
                        color='-log10P',
                        size_max=40,
                        color_continuous_scale='RdYlBu_r',
                        title='KEGG Pathway Enrichment Bubble Chart (Top 30 by P-value)',
                        labels={
                            'Rich_Factor': 'Rich Factor',
                            'label': 'Pathway (KEGG ID: Name)',
                            'Candidate_Hits': '# Compounds',
                            '-log10P': '-log10(P-value)',
                        },
                        hover_data={'KEGG_ID': True, 'Candidate_Hits': True, 'P_value': ':.2e'},
                    )
                    fig_bubble.update_layout(height=max(400, len(plot_df) * 18))
                    fig_bubble.update_yaxes(tickfont=dict(size=9))
                    st.plotly_chart(fig_bubble, use_container_width=True)

                    # ---- 柱状图（按富集数排序）----
                    st.markdown("**Top Pathways by Hit Count**")
                    bar_df = sig_df.nlargest(20, 'Candidate_Hits')
                    bar_df['label'] = bar_df['KEGG_ID'] + ': ' + bar_df['Pathway_Name']
                    bar_df['label'] = bar_df['label'].apply(
                        lambda x: x[:55] + '...' if len(str(x)) > 55 else x
                    )
                    fig_bar = px.bar(
                        bar_df,
                        x='Candidate_Hits',
                        y='label',
                        orientation='h',
                        color='-log10P',
                        color_continuous_scale='Viridis',
                        title='Top 20 Enriched Pathways by Hit Count',
                        labels={
                            'Candidate_Hits': '# Compounds in Pathway',
                            'label': 'Pathway',
                            '-log10P': '-log10(P-value)',
                        },
                        hover_data={'Rich_Factor': True, 'P_value': ':.2e'},
                    )
                    fig_bar.update_layout(height=max(400, len(bar_df) * 18), yaxis={'autorange': 'reversed'})
                    fig_bar.update_yaxes(tickfont=dict(size=9))
                    st.plotly_chart(fig_bar, use_container_width=True)

                    # ---- 结果表格 ----
                    st.markdown("**Full Enrichment Results (sorted by P-value)**")
                    display_cols = ['KEGG_ID', 'Pathway_Name', 'Candidate_Hits', 'Rich_Factor',
                                    'P_value', 'P_value_BH', 'Hit_Metabolites']
                    st.dataframe(
                        sig_df[display_cols].style.format({
                            'Rich_Factor': '{:.3f}',
                            'P_value': '{:.2e}',
                            'P_value_BH': '{:.2e}',
                        }, na_rep='-'),
                        use_container_width=True,
                        height=400,
                    )

                    # ---- KEGG Mapper 链接 ----
                    st.markdown("**KEGG Mapper Links**")
                    st.markdown(
                        "点击下方链接，在 KEGG Mapper 中将自动高亮显示命中的化合物："
                    )
                    unique_pids = sig_df['KEGG_ID'].dropna().unique()
                    top_pids = list(unique_pids[:10])

                    # 方案1：KEGG Mapper search_pathway 模式（以通路ID为关键词搜索）
                    mapper_url = (
                        "https://www.genome.jp/kegg/mapper/search_pathway/"
                        + "?keywords=" + '+'.join([str(p) for p in top_pids])
                    )

                    # 方案2：KEGG Search（搜索化合物）
                    kegg_search_url = (
                        "https://www.genome.jp/kegg/compound/"
                        + "?search=" + '+'.join([str(p) for p in top_pids[:5]])
                    )

                    # 方案3：KEGG BRITE（浏览生物通路）
                    brite_url = "https://www.genome.jp/kegg/brite.html"

                    st.markdown("**KEGG Pathway Visualization Links**")
                    st.markdown(
                        f"[1. KEGG Mapper - Search Pathways (Top 10 IDs)]({mapper_url})",
                        unsafe_allow_html=True
                    )
                    st.markdown(
                        f"[2. KEGG Search - Search by Pathway IDs]({kegg_search_url})",
                        unsafe_allow_html=True
                    )
                    st.markdown(
                        f"[3. KEGG BRITE - Browse Biological Pathways]({brite_url})",
                        unsafe_allow_html=True
                    )
                    st.caption(
                        "提示：KEGG Mapper 需手动在高亮页面输入化合物名称。"
                        "MetaboAnalyst（https://metaboanalyst.ca）支持直接上传化合物列表进行通路富集可视化。"
                    )

                    # ---- 下载 ----
                    bytes_enr = df_to_excel_bytes(sig_df[display_cols])
                    st.download_button(
                        "Download Enrichment Results (P<0.05)",
                        data=bytes_enr,
                        file_name="pathway_enrichment_results.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                    )
                else:
                    st.warning(
                        "候选化合物中未找到带 KEGG pathway_id 的注释。"
                        "请检查注释文件是否包含有效的 pathway_id。"
                    )

            # ---- 通路分布概览（不依赖按钮）----
            st.markdown("---")
            st.subheader("Pathway Coverage Overview")
            if st.session_state.get('pharma_df') is not None:
                pharma_df = st.session_state['pharma_df']
                # 提取所有 pathway_id 并统计
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

                    fig_dist = px.bar(
                        pw_display,
                        x='KEGG_ID',
                        y='Count',
                        color='Count',
                        color_continuous_scale='Blues',
                        title='KEGG Pathway ID Distribution (Top 30)',
                        labels={'KEGG_ID': 'KEGG Pathway ID', 'Count': 'Compound Count'},
                        text='Pathway_Name',
                    )
                    fig_dist.update_layout(height=450, xaxis={'tickangle': 45})
                    st.plotly_chart(fig_dist, use_container_width=True)
                else:
                    st.info("No KEGG pathway IDs found in pharmacological candidates.")

        # ---------- Tab 5: 药理数据库查询 ----------
        with tab5:
            st.subheader("Pharmacology Database Integration")
            st.markdown(
                "上传 TCMSP / DrugBank / TTD 离线数据库文件，"
                "为候选化合物匹配药理活性记录，并计算**药理证据评分**（纳入明星分子评分体系）。"
            )

            if match_pharma_db is None:
                st.error(
                    "pharma_db.py 模块未正确加载。"
                    "请确保 pharma_db.py 与 app.py 在同一目录下。"
                )
            else:
                # 文件上传区
                db_col1, db_col2, db_col3 = st.columns(3)
                with db_col1:
                    tcmsp_file = st.file_uploader(
                        "上传 TCMSP 数据库文件",
                        type=['xlsx', 'xls', 'csv', 'tsv'],
                        help="TCMSP (Traditional Chinese Medicine Systems Pharmacology) 数据库，"
                             "包含 OB(%)、DL、BBB 等药代动力学参数。"
                    )
                    if tcmsp_file:
                        st.success(f"Loaded: {tcmsp_file.name}")
                with db_col2:
                    drugbank_file = st.file_uploader(
                        "上传 DrugBank 数据库文件",
                        type=['csv', 'tsv', 'xml'],
                        help="DrugBank 数据库，包含药物靶点、适应症、药理活性信息。"
                    )
                    if drugbank_file:
                        st.success(f"Loaded: {drugbank_file.name}")
                with db_col3:
                    ttd_file = st.file_uploader(
                        "上传 TTD 数据库文件",
                        type=['xlsx', 'xls', 'csv', 'tsv'],
                        help="TTD (Therapeutic Target Database)，包含药物靶点与疾病关联信息。"
                    )
                    if ttd_file:
                        st.success(f"Loaded: {ttd_file.name}")

                use_stitch = st.checkbox(
                    "查询 STITCH 化学相互作用数据库（需要网络连接）",
                    value=False,
                    help="STITCH (http://stitch.embl.de) 提供化合物-蛋白质相互作用信息，"
                         "查询需要能够访问国际互联网。"
                )

                if st.button("Run Pharma DB Query", type="primary", use_container_width=True):
                    if load_tcmsp is None:
                        st.error("pharma_db 模块未安装")
                    else:
                        with st.spinner("Loading database files..."):
                            tcmsp_data = load_tcmsp(tcmsp_file) if tcmsp_file else {}
                            drugbank_data = load_drugbank(drugbank_file) if drugbank_file else {}
                            ttd_data = load_ttd(ttd_file) if ttd_file else {}

                            st.info(
                                f"TCMSP: {len(tcmsp_data)} compounds | "
                                f"DrugBank: {len(drugbank_data)} entries | "
                                f"TTD: {len(ttd_data)} entries"
                            )

                        with st.spinner("Matching compounds to pharmacology databases..."):
                            metabolites = st.session_state['pharma_df']['Metabolite'].tolist()

                            def progress_callback(current, total):
                                progress_bar.progress(int(current / total * 100))

                            progress_bar = st.progress(0)
                            pharma_match_df = match_pharma_db(
                                metabolites=metabolites,
                                tcmsp_data=tcmsp_data,
                                drugbank_data=drugbank_data,
                                ttd_data=ttd_data,
                                use_stitch=use_stitch,
                                progress_callback=progress_callback,
                            )
                            progress_bar.empty()
                            st.session_state['pharma_match_df'] = pharma_match_df

                            # 统计命中情况
                            n_tcmsp = pharma_match_df['TCMSP_OB'].notna().sum()
                            n_drugbank = (pharma_match_df['DrugBank_Targets'] != '-').sum()
                            n_ttd = (pharma_match_df['TTD_Targets'] != '-').sum()
                            n_active = (pharma_match_df['Pharma_Evidence_Score'] > 0).sum()
                            st.success(
                                f"匹配完成！TCMSP命中: {n_tcmsp} | DrugBank命中: {n_drugbank} | "
                                f"TTD命中: {n_ttd} | 有药理证据: {n_active} / {len(pharma_match_df)}"
                            )

                        # 保存匹配结果用于更新明星分子评分
                        st.session_state['pharma_match_df'] = pharma_match_df

                        # 更新明星分子评分（加入药理证据维度）
                        if 'star_df' in st.session_state and len(st.session_state['star_df']) > 0:
                            with st.spinner("Recalculating star scores with pharma evidence..."):
                                star_df = st.session_state['star_df'].copy()
                                # 合并药理证据评分
                                score_map = pharma_match_df.set_index('Metabolite')['Pharma_Evidence_Score'].to_dict()
                                star_df['_pharma_evidence'] = star_df['Metabolite'].map(score_map).fillna(0)

                                # 更新综合评分：原三维度 + 药理证据
                                # 重新归一化药理证据分数
                                pe_max = star_df['_pharma_evidence'].max()
                                star_df['_pharma_evidence_norm'] = (
                                    star_df['_pharma_evidence'] / pe_max if pe_max > 0
                                    else star_df['_pharma_evidence'] * 0
                                )

                                # 调整权重：VIP×0.25 + diff×0.25 + pharma_class×0.30 + pharma_evidence×0.20
                                star_df['_star_score_v2'] = (
                                    star_df['_vip_norm'].fillna(0) * 0.25 +
                                    star_df['_diff_norm'].fillna(0) * 0.25 +
                                    star_df['_pharma_score'].fillna(0) * 0.30 +
                                    star_df['_pharma_evidence_norm'].fillna(0) * 0.20
                                )
                                star_df = star_df.sort_values('_star_score_v2', ascending=False).reset_index(drop=True)
                                st.session_state['star_df_v2'] = star_df
                            st.success("明星分子评分已更新（含药理证据维度）！请切换到 Star Molecules 标签页查看。")

                # 显示药理匹配结果
                if 'pharma_match_df' in st.session_state and len(st.session_state['pharma_match_df']) > 0:
                    pm_df = st.session_state['pharma_match_df']

                    st.markdown("---")
                    st.subheader("Pharmacology Match Results")

                    # 统计卡片
                    pc1, pc2, pc3, pc4 = st.columns(4)
                    n_tcmsp = pm_df['TCMSP_OB'].notna().sum()
                    n_db = (pm_df['DrugBank_Targets'] != '-').sum()
                    n_ttd = (pm_df['TTD_Targets'] != '-').sum()
                    n_active = (pm_df['Pharma_Evidence_Score'] > 0).sum()
                    pc1.metric("TCMSP Hits", f"{n_tcmsp}")
                    pc2.metric("DrugBank Hits", f"{n_db}")
                    pc3.metric("TTD Hits", f"{n_ttd}")
                    pc4.metric("With Evidence", f"{n_active}/{len(pm_df)}")

                    # 药理证据评分分布图
                    st.markdown("**Pharma Evidence Score Distribution**")
                    fig_ev = px.histogram(
                        pm_df,
                        x='Pharma_Evidence_Score',
                        nbins=20,
                        title='Pharmacological Evidence Score Distribution',
                        labels={'Pharma_Evidence_Score': 'Evidence Score', 'count': 'Count'},
                        color_discrete_sequence=['#2ECC71'],
                    )
                    fig_ev.update_layout(height=350)
                    st.plotly_chart(fig_ev, use_container_width=True)

                    # 结果表格
                    st.markdown("**Full Match Results**")
                    display_cols = [
                        'Metabolite', 'TCMSP_OB', 'TCMSP_DL',
                        'DrugBank_Targets', 'DrugBank_Indications',
                        'TTD_Targets', 'TTD_Diseases',
                        'Pharma_Evidence_Score', 'Pharma_Evidence'
                    ]
                    avail_cols = [c for c in display_cols if c in pm_df.columns]
                    st.dataframe(
                        pm_df[avail_cols].style.format({
                            'TCMSP_OB': lambda x: f'{x:.1f}' if pd.notna(x) else '-',
                            'TCMSP_DL': lambda x: f'{x:.3f}' if pd.notna(x) else '-',
                            'Pharma_Evidence_Score': '{:.2f}',
                        }, na_rep='-'),
                        use_container_width=True,
                        height=500,
                    )

                    # 下载匹配结果
                    match_bytes = io.BytesIO()
                    with pd.ExcelWriter(match_bytes, engine='openpyxl') as writer:
                        pm_df.to_excel(writer, index=False, sheet_name='Pharma_Match')
                    match_bytes.seek(0)
                    st.download_button(
                        "Download Pharma Match Results",
                        data=match_bytes.getvalue(),
                        file_name="pharma_database_match.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                    )
    else:
        # 未运行分析时显示说明
        st.info("Please upload files in the sidebar and click **Start Analysis** to begin.")
        st.markdown("---")
        st.subheader("Workflow")
        st.markdown("""
        **Step 1: Upload Files**
        - Upload all `*_diff.exp.xls` comparison files (multiple files supported)
        - Upload the `anno.xls` annotation file

        **Step 2: Configure Parameters**
        - Adjust VIP and P-value thresholds in the sidebar

        **Step 3: Analysis**
        - Click **Start Analysis** to run the pipeline

        **Step 4: Results**
        - 📊 **Data Overview**: Download Excel files and preview data
        - 📈 **Visualization**: Clustered heatmap + boxplots + VIP distribution
        - ⭐ **Star Molecules**: Filter and rank pharmacological candidates
        - 🧬 **KEGG Enrichment**: Pathway enrichment analysis and bubble chart
        """)


if __name__ == '__main__':
    main()
