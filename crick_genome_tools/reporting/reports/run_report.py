import argparse
import os
import sys

from crick_genome_tools.reporting.reports.vector_core_aav_report import VectorCoreAavReport

def run(data_path, report_type):
    print(f"Running with data_path: {data_path}, report_type: {report_type}")

    report = VectorCoreAavReport(data_path)
    report.generate_report([])

if __name__ == "__main__":
    # Check args
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--report_type", required=True)
    args, _ = parser.parse_known_args()
    data_path = args.data_path
    report_type = args.report_type

    # Check data path exists
    if not os.path.exists(data_path):
        print(f"ERROR: Data path does not exist: {data_path}")
        sys.exit(1)

    # Run
    run(data_path, report_type)
