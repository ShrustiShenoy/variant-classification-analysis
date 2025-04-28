# variant-classification-analysis
This repository contains a Python script (variant_classification_all_stages3.py) for classifying genomic variants as coding or non-coding across multiple stages of esophageal squamous cell carcinoma (ESCC). The script processes variant data from an Excel file, maps variants to genomic regions using a GTF annotation file, and validates classifications against provided annotations. It generates detailed Excel outputs and a bar graph visualization.

# Overview

The variant_classification_all_stages3.py script analyzes single nucleotide polymorphisms (SNPs) from ESCC samples across four stages (I–IV). It uses the pybedtools library to map variants to coding (CDS, exons) and non-coding (UTRs) regions based on a GRCh38 GTF annotation file. The script classifies variants as "Coding," "Non-Coding," or "Unknown," validates these against the input Variant_Classification field, and logs discrepancies. Results are saved in an Excel file with per-stage summaries and a global summary, accompanied by a bar graph visualization.This script is designed for researchers and bioinformaticians studying cancer genomics, particularly those working with ESCC variant data.

# Features
Variant Mapping: Maps variants to CDS, exon, and UTR regions using a GRCh38 GTF file.

Hierarchical Classification: Classifies variants as Coding (CDS/exon overlap) or Non-Coding (UTR or no overlap), with optional transcript ID matching.

Validation: Compares computed classifications with expected classifications from Variant_Classification

# Output
Excel file with per-stage variant data, summaries, global summary, and discrepancy logs.
Bar graph comparing coding, non-coding, and unknown variants across stages.

# Prerequisites
# Dependencies
pip install pandas pybedtools openpyxl matplotlib
# system requirements
sudo apt update
sudo apt install bedtools
# Usuage
python3 variant_classification_all_stages3.py

