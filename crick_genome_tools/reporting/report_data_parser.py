"""
Class for parsing report data.
"""

import json
import logging
import os
import pickle
from enum import Enum

import numpy as np
import pandas as pd

from crick_genome_tools.io.vcf import generate_merged_vcf_report
from crick_genome_tools.reporting.custom.mosdepth_parser import parse_mosdepth_per_base
from crick_genome_tools.reporting.custom.samtools_parser import parse_samtools_flagstat, parse_samtools_idxstats
from crick_genome_tools.reporting.tqc.configuration import ToulligqcConf
from crick_genome_tools.reporting.tqc.fastq_extractor import FastqExtractor


log = logging.getLogger(__name__)


class ReportDataParser:
    def __init__(self, data_folder):
        self.data_folder = data_folder
        self.result_dict = {}
        self.dataframe_dict = {}
        self.merged_dataframe_dict = {}
        self.summary_data = None

        # List dir to get folders only
        self.folder_names = os.listdir(data_folder)
        self.folder_names = [folder_name for folder_name in self.folder_names if os.path.isdir(os.path.join(data_folder, folder_name))]

    def get_data(self, vcf_tools=[]):
        """
        Get data from all sources.
        """
        # Get summary data
        summary_file = os.path.join(self.data_folder, "summary.json")
        if os.path.exists(summary_file):
            with open(summary_file, "r", encoding="UTF-8") as f:
                self.summary_data = json.load(f)

        # Process each folder
        for folder_name in self.folder_names:
            folder_path = os.path.join(self.data_folder, folder_name)
            log.info(f"Processing folder: {folder_name}")

            # Switch data source
            if folder_name == "toulligqc":
                log.info("Processing toulligqc data")
                self.get_toulligqc_data(folder_path)
            elif folder_name == "samtools_host":
                log.info("Processing samtools host data")
                self.get_samtools_flagstat_data(folder_path, ".host", "host")
            elif folder_name == "samtools_contaminent":
                log.info("Processing samtools contaminent data")
                self.get_samtools_contam_data(folder_path, ".contam", "contam")
            elif folder_name == "samtools_alignment":
                log.info("Processing samtools alignment data")
                self.get_samtools_flagstat_data(folder_path, ".viral", "align")
            elif folder_name == "coverage":
                log.info("Processing coverage data")
                self.get_mosdepth_data(folder_path)
            elif folder_name == "ref":
                log.info("Processing reference data")
                self.get_ref_data(folder_path)
            elif folder_name == "consensus":
                log.info("Processing consensus data")
                self.get_consensus_data(folder_path)
            elif folder_name == "variants":
                log.info("Processing variant data")
                self.get_variant_data(folder_path, vcf_tools)
            elif folder_name == "variants_compressed":
                log.info("Processing compressed variant data")
                self.get_compressed_variant_data(folder_path)
            elif folder_name == "annotation":
                log.info("Processing annotation data")
                self.get_annotation_data(folder_path)
            elif folder_name == "samplesheet":
                log.info("Processing samplesheet data")
                self.get_samplesheet_data(folder_path)
            elif folder_name == "count_table":
                log.info("Processing count table data")
                self.get_count_table_data(folder_path)
            elif folder_name == "truncation":
                log.info("Processing truncation data")
                self.get_truncation_data(folder_path)
            else:
                log.error(f"Unknown folder: {folder_name}")

        # Sort all dictionaries by sample_id
        self.result_dict = dict(sorted(self.result_dict.items()))
        self.dataframe_dict = dict(sorted(self.dataframe_dict.items()))

        log.info("Report parsing complete")

    def save_data(self, output_file):
        """
        Save data to pickle file.
        """
        with open(output_file, "wb") as f:
            pickle.dump(self, f)

    def get_toulligqc_data(self, folder_path):
        """
        Get data from toulligqc.
        """

        # Default config
        config = ToulligqcConf()
        config["threshold"] = "10"
        config["batch_size"] = "1000"
        config["thread"] = "10"
        config["barcoding"] = "False"
        config["quiet"] = "False"
        config["barcode_selection"] = []

        # Get fastq files
        fastq_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".fastq.gz")]
        if len(fastq_files) == 0:
            log.error(f"No fastq files found in {folder_path}")
            return
        # Process each fastq file
        for fastq_file in fastq_files:
            config["fastq"] = os.path.join(folder_path, fastq_file)

            # Check for sample_id entry
            sample_id = fastq_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            if sample_id not in self.dataframe_dict:
                self.dataframe_dict[sample_id] = {}

            # Extract data
            extractor = FastqExtractor(config)
            extractor.init()
            result_dict = {}
            extractor.extract(result_dict)

            # Extract sample id from the first part of the fastq file name
            self.result_dict[sample_id]["toulligqc"] = result_dict
            self.dataframe_dict[sample_id]["toulligqc"] = extractor.dataframe_dict
            log.info(f"Processed fastq file: {fastq_file}")

    def get_samtools_flagstat_data(self, folder_path, clean_ext, data_suffix):
        """
        Get data from samtools reports
        """
        # Parse data
        data_dict = parse_samtools_flagstat(folder_path, clean_ext)

        # Extract summary table for primary reads
        rows = []
        for sample, metrics in data_dict.items():
            mapped = metrics["primary_mapped"]
            total = metrics["primary"]
            unmapped = total - mapped
            percent_mapped = round((mapped / total) * 100, 2)
            rows.append({"Sample": sample, "Total": total, "Mapped": mapped, "Unmapped": unmapped, "Percent mapped (%)": percent_mapped})

        df = pd.DataFrame(rows)
        df = df.sort_values(by=["Sample"])

        # Add data to merged
        self.merged_dataframe_dict["samtools_" + data_suffix] = df

    def get_samtools_contam_data(self, folder_path, clean_ext, data_suffix):
        """
        Get data from samtools reports
        """
        # Parse data
        data_dict = parse_samtools_idxstats(folder_path, clean_ext)

        # Extract summary table for primary reads
        rows = []
        for sample, metrics in data_dict.items():
            contigs = metrics.keys()
            row_data = {}
            row_data["Sample"] = sample
            for idx, contig in enumerate(contigs):
                column_name = contig
                if contig == "*":
                    column_name = "Unmapped"
                    row_data[column_name] = metrics[contig]["mapped"]
                elif idx == 0:
                    column_name = "Primary"
                    row_data[column_name] = metrics[contig]["reads"]
                else:
                    row_data[column_name] = metrics[contig]["reads"]
            rows.append(row_data)

        df = pd.DataFrame(rows)
        df = df.sort_values(by=["Sample"])

        # Add data to merged
        self.merged_dataframe_dict["samtools_" + data_suffix] = df

    def get_mosdepth_data(self, folder_path):
        """
        Get data from mosdepth reports.
        """
        self.dataframe_dict["coverage_per_base"] = parse_mosdepth_per_base(folder_path)

    def get_ref_data(self, folder_path):
        """
        Get data from reference files.
        """
        # Get reference files
        ref_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".fasta")]
        for ref_file in ref_files:
            sample_id = ref_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            # Read each line of the fasta file into a list
            with open(os.path.join(folder_path, ref_file), "r", encoding="UTF-8") as f:
                self.result_dict[sample_id]["ref"] = f.readlines()
            log.info(f"Processed reference file: {sample_id} - {ref_file}")

        # Get index files
        ref_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".fasta.fai")]
        for ref_file in ref_files:
            sample_id = ref_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            with open(os.path.join(folder_path, ref_file), "r", encoding="UTF-8") as f:
                self.result_dict[sample_id]["fai"] = f.readlines()
            log.info(f"Processed reference index file: {sample_id} - {ref_file}")

    def get_variant_data(self, folder_path, vcf_tools):
        """
        Get data from variant files.
        """
        # Get text variant files
        var_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".vcf")]
        for var_file in var_files:
            sample_id = var_file.split(".")[0]
            tool_name = var_file.split(".")[1]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            if "variants" not in self.result_dict[sample_id]:
                self.result_dict[sample_id]["variants"] = {}
            with open(os.path.join(folder_path, var_file), "r", encoding="UTF-8") as f:
                self.result_dict[sample_id]["variants"][tool_name] = f.readlines()
            log.info(f"Processed variant file: {sample_id} - {var_file}")

        if len(vcf_tools) > 0:
            # Make list of var files for each sample_id
            var_files_by_sample = {}
            for var_file in var_files:
                sample_id = var_file.split(".")[0]
                tool_name = var_file.split(".")[1]
                if tool_name in vcf_tools:
                    if sample_id not in var_files_by_sample:
                        var_files_by_sample[sample_id] = []
                    var_files_by_sample[sample_id].append(os.path.join(folder_path, var_file))

            # Process each sample_id
            for sample_id in var_files_by_sample.keys():
                if sample_id not in self.dataframe_dict:
                    self.dataframe_dict[sample_id] = {}
                var_files_by_sample[sample_id].sort(key=lambda x: vcf_tools.index(x.split(".")[1]))
                variants, header, processed_variants = generate_merged_vcf_report(var_files_by_sample[sample_id], vcf_tools)
                self.dataframe_dict[sample_id]["variants"] = pd.DataFrame(processed_variants, columns=header)
                log.info(f"Generated merged vcf report: {sample_id} - {var_files_by_sample[sample_id]}")

    def get_compressed_variant_data(self, folder_path):
        # Get binary variant files
        var_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".vcf.gz")]
        for var_file in var_files:
            sample_id = var_file.split(".")[0]
            tool_name = var_file.split(".")[1]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            if "variants_gz" not in self.result_dict[sample_id]:
                self.result_dict[sample_id]["variants_gz"] = {}
            with open(os.path.join(folder_path, var_file), "rb") as f:
                self.result_dict[sample_id]["variants_gz"][tool_name] = f.read()
            log.info(f"Processed variant gzip file: {sample_id} - {var_file}")

        # Get Tabix files
        var_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".vcf.gz.tbi")]
        for var_file in var_files:
            sample_id = var_file.split(".")[0]
            tool_name = var_file.split(".")[1]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            if "variants_tbi" not in self.result_dict[sample_id]:
                self.result_dict[sample_id]["variants_tbi"] = {}
            with open(os.path.join(folder_path, var_file), "rb") as f:
                self.result_dict[sample_id]["variants_tbi"][tool_name] = f.read()
            log.info(f"Processed variant tabix file: {sample_id} - {var_file}")

    def get_annotation_data(self, folder_path):
        # Get annotation files
        ann_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".gff")]
        for ann_file in ann_files:
            sample_id = ann_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            # Read each line of the fasta file into a list
            with open(os.path.join(folder_path, ann_file), "r", encoding="UTF-8") as f:
                self.result_dict[sample_id]["annotation"] = f.readlines()
            log.info(f"Processed annotation file: {sample_id} - {ann_file}")

    def get_samplesheet_data(self, folder_path):
        # Get path of first csv file in folder
        csv_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".csv")]
        self.merged_dataframe_dict["samplesheet"] = pd.read_csv(os.path.join(folder_path, csv_files[0]), encoding="UTF-8")
        log.info(f"Processed samplesheet file: {csv_files[0]}")

    def get_count_table_data(self, folder_path):
        csv_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".csv")]
        for csv_file in csv_files:
            sample_id = csv_file.split(".")[0]
            if sample_id not in self.dataframe_dict:
                self.dataframe_dict[sample_id] = {}
            self.dataframe_dict[sample_id]["count_table"] = pd.read_csv(os.path.join(folder_path, csv_file), encoding="UTF-8", sep="\t")
            log.info(f"Processed count table file: {csv_file}")

    def get_consensus_data(self, folder_path):
        # Get consensus files
        cons_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".fasta")]
        for cons_file in cons_files:
            sample_id = cons_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            # Read each line of the fasta file into a list
            with open(os.path.join(folder_path, cons_file), "r", encoding="UTF-8") as f:
                self.result_dict[sample_id]["consensus"] = f.readlines()
            log.info(f"Processed consensus file: {sample_id} - {cons_file}")

    def get_truncation_data(self, folder_path):
        bam_info_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".tsv")]
        for bam_info_file in bam_info_files:
            sample_id = bam_info_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            # Read the bam info and save the starts/ends positions
            bam_info_df = pd.read_csv(
                os.path.join(folder_path, bam_info_file),
                sep='\t',
                usecols=['Pos', 'EndPos'],
                dtype={
                    'Pos': np.uint32,
                    'EndPos': np.uint32
                })
            bam_info_df = bam_info_df.rename(columns={'Pos': 'Read Start', 'EndPos': 'Read End'})
            self.result_dict[sample_id]["truncation"] = bam_info_df
            log.info(f"Processed truncation file: {sample_id} - {bam_info_file}")

            bam_info_df = pd.read_csv(
                os.path.join(folder_path, bam_info_file),
                sep='\t',
                usecols=['Ref', 'Read', 'Pos', 'EndPos', 'ReadLen', 'Strand', 'IsSec', 'IsSup'],
                dtype = {
                'Read': 'string',
                'Ref': 'string',
                'Pos': np.int32,
                'EndPos': np.int32,
                'ReadLen': np.int32,
                'Strand': np.int8,
                'IsSec': np.int8,
                'IsSup': np.int8,
                })
            bam_info_df[['Strand', 'IsSec', 'IsSup']] = bam_info_df[['Strand', 'IsSec', 'IsSup']].astype(np.bool_)

            itr_length = 130
            itr_fl_threshold = 20
            payload_threshold = 100
            itr1_starts = 0
            itr1_ends = itr_length
            itr2_ends = bam_info_df['EndPos'].max()
            itr2_starts = itr2_ends - itr_length

            starts = bam_info_df['Pos']
            ends = bam_info_df['EndPos']
            itr1_full = itr1_starts + itr_fl_threshold
            itr2_full = itr2_ends - itr_fl_threshold
            payload5_full = itr1_ends + payload_threshold
            payload3_full = itr2_starts - payload_threshold

            full_5prime = (starts >= itr1_starts) & (starts < itr1_full)
            full_3prime = (ends > itr2_full) & (ends <= itr2_ends)
            partial_5prime = (starts >= itr1_full) & (starts <= itr1_ends)
            partial_3prime = (ends >= itr2_starts) & (ends <= itr2_full)
            full_payload = (starts <= payload5_full) & (ends >= payload3_full)
            starts_in_midsection = (starts > itr1_ends) & (starts <= itr2_starts)
            ends_in_midsection = (ends > itr1_ends) & (ends < itr2_starts)

            conditions = [
                full_5prime & full_3prime,
                partial_5prime & full_3prime,
                partial_5prime & partial_3prime,
                full_5prime & partial_3prime,
                full_payload,
                full_5prime & ends_in_midsection,
                starts_in_midsection & full_3prime,
                partial_5prime & ends_in_midsection,
                starts_in_midsection & partial_3prime,
                (starts >= itr1_starts) & (ends <= itr1_ends),
                (starts >= itr2_starts) & (ends <= itr2_ends),
                (starts < itr1_starts) & (ends > itr2_ends),
                (starts < itr1_starts) & (ends >= itr1_starts),
                (starts <= itr2_ends) & (ends > itr2_ends),
                (starts < itr1_starts) & (ends < itr1_starts),
                (starts > itr2_ends) & (ends > itr2_ends)
            ]

            choices = [
                AlnType.complete,
                AlnType.par5_full3,
                AlnType.par5_par3,
                AlnType.full5_par3,
                AlnType.full_payload,
                AlnType.full5_par_mid,
                AlnType.par_mid_full3,
                AlnType.par5_par_mid,
                AlnType.par_mid_par3,
                AlnType.itr5_only,
                AlnType.itr3_only,
                AlnType.ext_itr,
                AlnType.vec_bb_5,
                AlnType.vec_bb_3,
                AlnType.bb,
                AlnType.bb
            ]

            choices_simple = [
                AlnType.complete,
                AlnType.full_payload,
                AlnType.full_payload,
                AlnType.full_payload,
                AlnType.full_payload,
                AlnType.truncated_payload,
                AlnType.truncated_payload,
                AlnType.truncated_payload,
                AlnType.truncated_payload,
                AlnType.itr5_only,
                AlnType.itr3_only,
                AlnType.ext_itr,
                AlnType.vec_bb_5,
                AlnType.vec_bb_3,
                AlnType.bb,
                AlnType.bb
            ]

            bam_info_df["aln_type"] = np.select(
                conditions,
                choices,
                default=AlnType.unknown
            )
            self.result_dict[sample_id]["truncation_type"] = bam_info_df

            bam_info_df_simple = bam_info_df.copy()
            bam_info_df_simple["aln_type"] = np.select(
                conditions,
                choices_simple,
                default=AlnType.unknown
            )
            self.result_dict[sample_id]["truncation_type_simple"] = bam_info_df_simple

class AlnType(str, Enum):
    """Enum for Assigning categories to alignments.

    An alignment category defines its ITR and midsection status as well as whether
    the alignment maps to the vector backbone.

    Subclassing str allows us to access the values as strings and not have to
    do .value all over the place.
    """

    # These are alignments that represent almost full AAV genomes. They have varying
    # amounts of ITR on both sides of the alignment and contain full mid-sections
    complete = 'Complete'
    full5_par3 = 'Full 5` ITR and partial 3` ITR'
    par5_full3 = 'Partial 5` ITR and full 3` ITR'
    par5_par3 = 'Partial 5` ITR and partial 3` ITR'
    full_payload = 'Full payload'

    # These alignments are truncated at the mid-section region but contain some
    # ITR region on one of the ends
    truncated_payload = 'Truncated payload'
    full5_par_mid = 'Full 5` ITR and partial payload'
    par_mid_full3 = 'Partial payload section and full 3` ITR'
    par5_par_mid = 'Partial 5` ITR and partial payload'
    par_mid_par3 = 'Partial payload and partial 3` ITR'

    # Alignment starts and ends within ITR
    itr5_only = '5` ITR'
    itr3_only = '3` ITR'

    # Transgene plamsid backbone alignments
    vec_bb_5 = 'Vector backbone - 5` ends'
    vec_bb_3 = 'Vector backbone - 3` ends'
    bb = 'Backbone'

    ext_itr = 'Extended ITR-ITR region'
    unknown = 'Unknown'

    def __str__(self):
        return self.value
