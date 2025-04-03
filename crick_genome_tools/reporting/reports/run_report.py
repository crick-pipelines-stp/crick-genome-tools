"""
Run a report on a dataset.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring,invalid-name

import argparse
import os
import sys
import uuid

from crick_genome_tools.reporting.reports.vector_core_aav_report import VectorCoreAavReport

STATIC_TMP_DIR = "crick_genome_tools/reporting/reports/static/tmp"


def run(data_path, report_type):
    if report_type == "aav":
        file_prefix = str(f"{STATIC_TMP_DIR}/{uuid.uuid4().hex}")
        report = VectorCoreAavReport("20241119_1033_P2S-02348-A_PAY61327_7eb45677", data_path, tmp_dir=file_prefix)
        report.generate_report()
    else:
        print(f"ERROR: Unknown report type: {report_type}")
        sys.exit(1)

if __name__ == "__main__":
    # Check args
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--report_type", required=True)
    args, _ = parser.parse_known_args()

    # Check data path exists
    if not os.path.exists(args.data_path):
        print(f"ERROR: Data path does not exist: {args.data_path}")
        sys.exit(1)

    # Run
    run(args.data_path, args.report_type)
