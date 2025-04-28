import pandas as pd
import pybedtools
from pybedtools import BedTool
import logging
import sys
from collections import defaultdict
import matplotlib.pyplot as plt

# Set up logging for debugging and error tracking
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Step 1: Load the GTF annotation file
logger.info("Loading GTF annotation file...")
try:
    gtf_file = "annotation/Homo_sapiens.GRCh38.113.gtf/Homo_sapiens.GRCh38.113.gtf"
    gtf = BedTool(gtf_file)
except FileNotFoundError:
    logger.error("GTF file not found. Please ensure the file exists in the 'annotation' directory.")
    sys.exit(1)
except Exception as e:
    logger.error(f"Error loading GTF file: {e}")
    sys.exit(1)

# Filter GTF for CDS (strictly coding regions)
logger.info("Filtering GTF for CDS regions...")
cds_gtf = gtf.filter(lambda x: x[2] == "CDS").saveas()

# Filter GTF for exons (to refine coding status, will exclude UTRs later)
logger.info("Filtering GTF for exon regions...")
exon_gtf = gtf.filter(lambda x: x[2] == "exon").saveas()

# Filter GTF for UTR regions (non-coding)
logger.info("Filtering GTF for UTR regions...")
utr_gtf = gtf.filter(lambda x: x[2] in ["5'UTR", "3'UTR"]).saveas()

# Step 2: Load the Excel file and process each sheet
logger.info("Loading Excel file...")
excel_file = "ESCC_UP_analysisReport_2025-04-24_05-04-01.xlsx"
try:
    xl = pd.ExcelFile(excel_file)
except FileNotFoundError:
    logger.error("Excel file not found. Please ensure the file exists in the working directory.")
    sys.exit(1)
except Exception as e:
    logger.error(f"Error loading Excel file: {e}")
    sys.exit(1)

# Define the stages to process
stages = [f"Stage{i}_SNP_Position_Data" for i in ["I", "II", "III", "IV"]]
logger.info(f"Sheets to process: {stages}")

# Initialize a dictionary to store summaries for each stage
stage_summaries = defaultdict(dict)
discrepancy_logs = []
successful_stages = []

# Process each sheet
output_file = "variant_coding_classification_all_stages.xlsx"
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    for stage in stages:
        if stage not in xl.sheet_names:
            logger.warning(f"Sheet {stage} not found in Excel file. Skipping...")
            continue

        logger.info(f"Processing sheet: {stage}...")
        try:
            variant_df = pd.read_excel(excel_file, sheet_name=stage)
        except Exception as e:
            logger.error(f"Error reading sheet {stage}: {e}")
            continue

        # Check if the DataFrame is empty
        if variant_df.empty:
            logger.warning(f"No variants found in sheet {stage}. Skipping...")
            continue

        # Ensure required columns exist
        required_columns = ["Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Transcript_ID"]
        missing_columns = [col for col in required_columns if col not in variant_df.columns]
        if missing_columns:
            logger.warning(f"Missing columns in sheet {stage}: {missing_columns}. Proceeding without missing columns...")
            if "Transcript_ID" in missing_columns:
                logger.warning("Transcript_ID column missing. Classification will rely on coordinates only.")

        # Data quality checks
        variant_df = variant_df[variant_df["Start_Position"].notna() & variant_df["End_Position"].notna()]
        variant_df = variant_df[variant_df["Start_Position"] >= 0]
        variant_df = variant_df[variant_df["End_Position"] >= variant_df["Start_Position"]]

        # Standardize chromosome format
        variant_df["Chromosome"] = variant_df["Chromosome"].astype(str).str.replace("chr", "")

        # Step 3: Prepare variant data as a BED file for mapping
        logger.info(f"Converting variant data from {stage} to BED format...")
        variant_bed_df = variant_df[["Chromosome", "Start_Position", "End_Position"]].copy()
        variant_bed_df["Start_Position"] = variant_bed_df["Start_Position"].astype(int)
        variant_bed_df["End_Position"] = variant_bed_df["End_Position"].astype(int)
        variant_bed = BedTool.from_dataframe(variant_bed_df)

        # Step 4: Map variants to CDS regions (strictly coding)
        logger.info(f"Mapping variants from {stage} to CDS regions...")
        cds_intersection = variant_bed.intersect(cds_gtf, wa=True, wb=True)

        # Step 5: Map variants to exon regions
        logger.info(f"Mapping variants from {stage} to exon regions...")
        exon_intersection = variant_bed.intersect(exon_gtf, wa=True, wb=True)

        # Step 6: Map variants to UTR regions (non-coding)
        logger.info(f"Mapping variants from {stage} to UTR regions...")
        utr_intersection = variant_bed.intersect(utr_gtf, wa=True, wb=True)

        # Step 7: Classify variants as coding or non-coding with Transcript_ID
        logger.info(f"Classifying variants in {stage} as coding or non-coding...")
        cds_variants = {}
        for interval in cds_intersection:
            chrom, start = interval[0], int(interval[1])
            try:
                attributes = interval[8]
                if "transcript_id" in attributes:
                    transcript_id = attributes.split("transcript_id ")[1].split(";")[0].strip('"')
                else:
                    transcript_id = "N/A"
                    logger.warning(f"No transcript_id found in CDS attributes: {attributes}")
            except (IndexError, KeyError) as e:
                transcript_id = "N/A"
                logger.warning(f"Error parsing transcript_id from CDS attributes: {interval[8]}, error: {e}")
            cds_variants[(chrom, start)] = transcript_id

        exon_variants = {}
        for interval in exon_intersection:
            chrom, start = interval[0], int(interval[1])
            try:
                attributes = interval[8]
                if "transcript_id" in attributes:
                    transcript_id = attributes.split("transcript_id ")[1].split(";")[0].strip('"')
                else:
                    transcript_id = "N/A"
                    logger.warning(f"No transcript_id found in exon attributes: {attributes}")
            except (IndexError, KeyError) as e:
                transcript_id = "N/A"
                logger.warning(f"Error parsing transcript_id from exon attributes: {interval[8]}, error: {e}")
            exon_variants[(chrom, start)] = transcript_id

        utr_variants = {}
        for interval in utr_intersection:
            chrom, start = interval[0], int(interval[1])
            try:
                attributes = interval[8]
                if "transcript_id" in attributes:
                    transcript_id = attributes.split("transcript_id ")[1].split(";")[0].strip('"')
                else:
                    transcript_id = "N/A"
                    logger.warning(f"No transcript_id found in UTR attributes: {attributes}")
            except (IndexError, KeyError) as e:
                transcript_id = "N/A"
                logger.warning(f"Error parsing transcript_id from UTR attributes: {interval[8]}, error: {e}")
            utr_variants[(chrom, start)] = transcript_id

        # Hierarchical classification with Transcript_ID
        def classify_variant(row):
            key = (str(row["Chromosome"]), int(row["Start_Position"]))
            transcript_id = row.get("Transcript_ID", None) if "Transcript_ID" in variant_df.columns else None

            if key in cds_variants:
                if transcript_id and cds_variants[key] == transcript_id:
                    return "Coding"  # CDS overlap with matching transcript
                elif not transcript_id:
                    return "Coding"  # CDS overlap without transcript check
            elif key in utr_variants:
                if transcript_id and utr_variants[key] == transcript_id:
                    return "Non-Coding"  # UTR overlap with matching transcript
                elif not transcript_id:
                    return "Non-Coding"  # UTR overlap without transcript check
            elif key in exon_variants:
                if transcript_id and exon_variants[key] == transcript_id:
                    return "Coding"  # Exon overlap with matching transcript
                elif not transcript_id:
                    return "Coding"  # Exon overlap without transcript check
            else:
                return "Non-Coding"  # No overlap with coding regions

            return "Unknown"  # Fallback for unmatched cases

        variant_df["Computed_Coding_Status"] = variant_df.apply(classify_variant, axis=1)

        # Step 8: Validate against Variant_Classification
        logger.info(f"Validating coding status against Variant_Classification in {stage}...")
        coding_classifications = {
            "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
            "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
            "Start_Codon", "Stop_Codon", "Nonstop_Mutation", "Splice_Region"
        }
        variant_df["Expected_Coding_Status"] = variant_df["Variant_Classification"].apply(
            lambda x: "Coding" if x in coding_classifications else "Non-Coding"
        )

        # Identify and log discrepancies
        variant_df["Discrepancy"] = variant_df["Computed_Coding_Status"] != variant_df["Expected_Coding_Status"]
        discrepant_rows = variant_df[variant_df["Discrepancy"]]
        for _, row in discrepant_rows.iterrows():
            key = (str(row["Chromosome"]), int(row["Start_Position"]))
            transcript_id = row.get("Transcript_ID", "N/A")
            discrepancy_logs.append({
                "Stage": stage,
                "Chromosome": row["Chromosome"],
                "Start_Position": row["Start_Position"],
                "End_Position": row["End_Position"],
                "Computed_Coding_Status": row["Computed_Coding_Status"],
                "Expected_Coding_Status": row["Expected_Coding_Status"],
                "Variant_Classification": row["Variant_Classification"],
                "Transcript_ID": transcript_id,
                "CDS_Transcript": cds_variants.get(key, "N/A"),
                "Exon_Transcript": exon_variants.get(key, "N/A"),
                "UTR_Transcript": utr_variants.get(key, "N/A")
            })

        # Step 9: Summarize the results for this stage
        logger.info(f"Summarizing results for {stage}...")
        computed_summary = variant_df["Computed_Coding_Status"].value_counts().to_dict()
        expected_summary = variant_df["Expected_Coding_Status"].value_counts().to_dict()
        discrepancies = variant_df["Discrepancy"].value_counts().to_dict()

        coding_count = computed_summary.get("Coding", 0)
        non_coding_count = computed_summary.get("Non-Coding", 0)
        unknown_count = computed_summary.get("Unknown", 0)
        discrepancy_count = discrepancies.get(True, 0)

        stage_summaries[stage] = {
            "Total_Variants": len(variant_df),
            "Computed_Coding": coding_count,
            "Computed_Non_Coding": non_coding_count,
            "Computed_Unknown": unknown_count,
            "Expected_Coding": expected_summary.get("Coding", 0),
            "Expected_Non_Coding": expected_summary.get("Non-Coding", 0),
            "Discrepancies": discrepancy_count
        }

        logger.info(f"Stage {stage} - Total Variants: {len(variant_df)}")
        logger.info(f"Stage {stage} - Computed Coding Variants: {coding_count}")
        logger.info(f"Stage {stage} - Computed Non-Coding Variants: {non_coding_count}")
        logger.info(f"Stage {stage} - Computed Unknown Variants: {unknown_count}")
        logger.info(f"Stage {stage} - Expected Coding Variants: {expected_summary.get('Coding', 0)}")
        logger.info(f"Stage {stage} - Expected Non-Coding Variants: {expected_summary.get('Non-Coding', 0)}")
        logger.info(f"Stage {stage} - Discrepancies: {discrepancy_count}")

        # Step 10: Save the results for this stage to the same Excel file
        logger.info(f"Saving results for {stage} to Excel file...")
        try:
            variant_df.to_excel(writer, sheet_name=stage, index=False)
            summary_df = pd.DataFrame({
                "Category": [
                    "Total Variants", "Computed Coding Variants", "Computed Non-Coding Variants",
                    "Computed Unknown Variants", "Expected Coding Variants", "Expected Non-Coding Variants",
                    "Discrepancies"
                ],
                "Count": [
                    len(variant_df), coding_count, non_coding_count, unknown_count,
                    expected_summary.get("Coding", 0), expected_summary.get("Non-Coding", 0),
                    discrepancy_count
                ]
            })
            summary_df.to_excel(writer, sheet_name=f"{stage}_Summary", index=False)
            successful_stages.append(stage)
        except Exception as e:
            logger.error(f"Error saving results for {stage} to Excel: {e}")

    # Step 11: Save a global summary across all stages
    if successful_stages:
        logger.info("Saving global summary across all stages...")
        global_summary_df = pd.DataFrame.from_dict(stage_summaries, orient="index")
        try:
            global_summary_df.to_excel(writer, sheet_name="Global_Summary")
            logger.info(f"Global summary saved to {output_file}")
        except Exception as e:
            logger.error(f"Error saving global summary to Excel: {e}")

        # Step 12: Save discrepancy logs
        logger.info("Saving discrepancy logs...")
        discrepancy_df = pd.DataFrame(discrepancy_logs)
        if not discrepancy_df.empty:
            try:
                discrepancy_df.to_excel(writer, sheet_name="Discrepancy_Log", index=False)
                logger.info("Discrepancy log saved to Discrepancy_Log sheet")
            except Exception as e:
                logger.error(f"Error saving discrepancy log to Excel: {e}")
    else:
        logger.error("No stages were processed successfully. No Excel output will be written.")
        # Create a dummy sheet to avoid the "At least one sheet must be visible" error
        pd.DataFrame({"Message": ["No data processed"]}).to_excel(writer, sheet_name="Error", index=False)

# Step 13: Create visualization (bar graph) if there are successful stages
if successful_stages:
    logger.info("Generating bar graph visualization...")
    stages_processed = [stage for stage in stages if stage in xl.sheet_names and stage in successful_stages]
    coding_counts = [stage_summaries[stage]["Computed_Coding"] for stage in stages_processed]
    non_coding_counts = [stage_summaries[stage]["Computed_Non_Coding"] for stage in stages_processed]
    unknown_counts = [stage_summaries[stage]["Computed_Unknown"] for stage in stages_processed]

    plt.figure(figsize=(10, 6))
    bar_width = 0.25
    index = range(len(stages_processed))

    bars1 = plt.bar(index, coding_counts, bar_width, label="Computed Coding Variants", color="skyblue")
    bars2 = plt.bar([i + bar_width for i in index], non_coding_counts, bar_width, label="Computed Non-Coding Variants", color="lightcoral")
    bars3 = plt.bar([i + 2 * bar_width for i in index], unknown_counts, bar_width, label="Computed Unknown Variants", color="lightgray")

    # Add numbers on top of bars
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2, height, f'{int(height)}', 
                     ha='center', va='bottom' if height > 0 else 'top')

    plt.xlabel("Stages")
    plt.ylabel("Count")
    plt.title("Computed Coding vs Non-Coding vs Unknown Variants by Stage")
    plt.xticks([i + bar_width for i in index], stages_processed)
    plt.legend()
    plt.tight_layout()
    plt.savefig("variant_coding_classification_bar_graph.png")
    logger.info("Bar graph saved as variant_coding_classification_bar_graph.png")

    # Print global summary to console
    print("\nGlobal Summary Across All Stages:")
    print(global_summary_df)
else:
    logger.warning("No successful stages to visualize.")