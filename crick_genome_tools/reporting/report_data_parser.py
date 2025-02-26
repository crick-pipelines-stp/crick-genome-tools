"""
Class for parsing report data.
"""

import logging
import os

from crick_genome_tools.reporting.toulligqc.configuration import ToulligqcConf
from crick_genome_tools.reporting.toulligqc.fastq_extractor import fastqExtractor

log = logging.getLogger(__name__)

class ReportDataParser:
    def __init__(self, data_folder):
        self.data_folder = data_folder
        self.result_dict = {}
        self.dataframe_dict = {}

        # List dir to get folders only
        self.folder_names = os.listdir(data_folder)
        self.folder_names = [folder_name for folder_name in self.folder_names if os.path.isdir(os.path.join(data_folder, folder_name))]

    def get_data(self):
        """
        Get data from all sources.
        """

        # Process each folder
        for folder_name in self.folder_names:
            folder_path = os.path.join(self.data_folder, folder_name)
            log.info(f"Processing folder: {folder_name}")

            # Switch data source
            if folder_name == "toulligqc":
                self.get_toulligqc_data(folder_path)
            else:
                log.error(f"Unknown folder: {folder_name}")

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
                self.dataframe_dict[sample_id] = {}

            # Extract data
            extractor = fastqExtractor(config)
            extractor.init()
            result_dict = {}
            extractor.extract(result_dict)

            # Extract sample id from the first part of the fastq file name
            self.result_dict[sample_id]["toulligqc"] = result_dict
            self.dataframe_dict[sample_id]["toulligqc"] = extractor.dataframe_dict
            log.info(f"Processed fastq file: {fastq_file}")
