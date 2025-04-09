"""
Class for parsing report data.
"""

import logging
import os
import pickle
import json

import pandas as pd

from crick_genome_tools.reporting.tqc.configuration import ToulligqcConf
from crick_genome_tools.reporting.tqc.fastq_extractor import FastqExtractor
from crick_genome_tools.reporting.custom.samtools_parser import parse_samtools_flagstat, parse_samtools_idxstats
from crick_genome_tools.reporting.custom.mosdepth_parser import parse_mosdepth_per_base
from crick_genome_tools.io.vcf import generate_merged_vcf_report


log = logging.getLogger(__name__)

class ReportDataParser:
    def __init__(self, data_folder):
        self.data_folder = data_folder
        self.result_dict = {}
        self.dataframe_dict = {}
        self.merged_dataframe_dict = {}
        self.summary_data = None

        # List dir to get folders only
        self.folder_names = os.listdir(data_folder)
        self.folder_names = [folder_name for folder_name in self.folder_names if os.path.isdir(os.path.join(data_folder, folder_name))]

    def get_data(self, vcf_tools = []):
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

            # Switch data source
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
            else:
                log.error(f"Unknown folder: {folder_name}")

        # Sort all dictionaries by sample_id
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

        # Default config
        config = ToulligqcConf()
        config["threshold"] = "10"
        config["batch_size"] = "1000"
        config["thread"] = "10"
        config["barcoding"] = "False"
        config["quiet"] = "False"
        config["barcode_selection"] = []

        # Get fastq files
        fastq_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".fastq.gz")]
        if len(fastq_files) == 0:
            log.error(f"No fastq files found in {folder_path}")
            return
        # Process each fastq file
        for fastq_file in fastq_files:
            config["fastq"] = os.path.join(folder_path, fastq_file)

            # Check for sample_id entry
            sample_id = fastq_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            if sample_id not in self.dataframe_dict:
                self.dataframe_dict[sample_id] = {}

            # Extract data
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
        # Parse data
        data_dict = parse_samtools_flagstat(folder_path, clean_ext)

        # Extract summary table for primary reads
        rows = []
        for sample, metrics in data_dict.items():
            mapped = metrics['primary_mapped']
            total = metrics['primary']
            unmapped = total - mapped
            percent_mapped = round((mapped / total) * 100, 2)
            rows.append({'Sample': sample, 'Total': total, 'Mapped': mapped, 'Unmapped': unmapped, 'Percent mapped (%)': percent_mapped})

        df = pd.DataFrame(rows)
        df = df.sort_values(by=['Sample'])

        # Add data to merged
        self.merged_dataframe_dict["samtools_" + data_suffix] = df

    def get_samtools_contam_data(self, folder_path, clean_ext, data_suffix):
        """
        Get data from samtools reports
        """
        # Parse data
        data_dict = parse_samtools_idxstats(folder_path, clean_ext)

        # Extract summary table for primary reads
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
        df = df.sort_values(by=['Sample'])

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
        # Get reference files
        ref_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".fasta")]
        for ref_file in ref_files:
            sample_id = ref_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            # Read each line of the fasta file into a list
            with open(os.path.join(folder_path, ref_file), "r", encoding="UTF-8") as f:
                self.result_dict[sample_id]["ref"] = f.readlines()
            log.info(f"Processed reference file: {sample_id} - {ref_file}")

        # Get index files
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
        # Get text variant files
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
            # Make list of var files for each sample_id
            var_files_by_sample = {}
            for var_file in var_files:
                sample_id = var_file.split(".")[0]
                tool_name = var_file.split(".")[1]
                if tool_name in vcf_tools:
                    if sample_id not in var_files_by_sample:
                        var_files_by_sample[sample_id] = []
                    var_files_by_sample[sample_id].append(os.path.join(folder_path, var_file))

            # Process each sample_id
            for sample_id in var_files_by_sample.keys():
                if sample_id not in self.dataframe_dict:
                    self.dataframe_dict[sample_id] = {}
                var_files_by_sample[sample_id].sort(key=lambda x: vcf_tools.index(x.split(".")[1]))
                variants, header, processed_variants = generate_merged_vcf_report(var_files_by_sample[sample_id], vcf_tools)
                self.dataframe_dict[sample_id]["variants"] = pd.DataFrame(processed_variants, columns=header)
                log.info(f"Generated merged vcf report: {sample_id} - {var_files_by_sample[sample_id]}")


    def get_compressed_variant_data(self, folder_path):
        # Get binary variant files
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
        # Get annotation files
        ann_files = [file_name for file_name in os.listdir(folder_path) if file_name.endswith(".gff")]
        for ann_file in ann_files:
            sample_id = ann_file.split(".")[0]
            if sample_id not in self.result_dict:
                self.result_dict[sample_id] = {}
            # Read each line of the fasta file into a list
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
            # Read each line of the fasta file into a list
            with open(os.path.join(folder_path, cons_file), "r", encoding="UTF-8") as f:
                self.result_dict[sample_id]["consensus"] = f.readlines()
            log.info(f"Processed consensus file: {sample_id} - {cons_file}")
