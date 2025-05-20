"""
Utilities for processing VCF files
"""

import logging

import pysam


log = logging.getLogger(__name__)


def determine_variant_type(ref, alt):
    """
    Determine the type of variant (SNP, INS, DEL, INDEL)

    Args:
        ref (str): Reference sequence
        alt (str): Alternate sequence

    Returns:
        str: Type of variant
    """
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) < len(alt) and ref == alt[: len(ref)]:
        return "INS"
    if len(ref) > len(alt) and alt == ref[: len(alt)]:
        return "DEL"
    return "INDEL"


def generate_merged_vcf_report(vcf_files: list, tool_names: list, output_file: str = None):
    """
    Generate a report from a list of VCF files. The first tool must be the consensus tool,
    the second tools must be the primary variant caller.

    Args:
        vcf_files (list): List of paths to VCF files
        tool_names (list): List of tool names corresponding to the VCF files
        output_file (str, optional): Path to the output file
    """
    # Load vcf files
    vcf_data = []
    for vcf_file in vcf_files:
        curr_vcf = pysam.VariantFile(vcf_file)  # pylint: disable=no-member
        vcf_data.append(curr_vcf)

    # Check we have at least 2 tools
    if len(vcf_data) < 2:
        raise ValueError("At least two VCF files are required for comparison")

    # Error if files and tool len mismatch
    if len(vcf_data) != len(tool_names):
        raise ValueError(f"Number of VCF files and tool names must match {tool_names} != {vcf_data}")

    # Loop vcf files and collect information
    variants = {}
    for idx, vcf_file in enumerate(vcf_data):
        tool = tool_names[idx]
        for record in vcf_file.fetch():
            # Extract basic information
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt_list = [str(a) for a in record.alts] if record.alts else []
            if len(alt_list) > 1:
                alt = ",".join(alt_list)
            elif len(alt_list) == 0:
                alt_list = ref
            var_type = determine_variant_type(ref, alt_list)
            qual = round(float(record.qual), 2)
            info = record.info

            # Init optional fields
            depth = None
            depth_fwd = None
            depth_rev = None
            depth_ref = None
            depth_alt = None
            depth_ref_fwd = None
            depth_ref_rev = None
            depth_alt_fwd = None
            depth_alt_rev = None
            allele_freq = None
            annotation = None

            # Tool specific info: medaka
            if tool == "medaka":
                depth = info["DP"]
                depth_fwd = info["DPS"][0]
                depth_rev = info["DPS"][1]

            # Tool specific info: claire3
            if tool == "clair3":
                sample = record.samples[0]
                depth = sample["DP"]
                depth_ref = sample["AD"][0]
                depth_alt = sample["AD"][1]
                allele_freq = sample["AF"][0]
                allele_freq = round(float(allele_freq), 4)

            # Tool specific info: lofreq
            if tool == "lofreq":
                depth = info["DP"]
                allele_freq = info["AF"]
                allele_freq = round(float(allele_freq), 4)
                depth_ref_fwd = info["DP4"][0]
                depth_ref_rev = info["DP4"][1]
                depth_alt_fwd = info["DP4"][2]
                depth_alt_rev = info["DP4"][3]

            # Tool specific info: freebayes
            if tool == "freebayes":
                depth = info["DP"]
                allele_freq = info["AF"]
                allele_freq = round(float(allele_freq[0]), 4)

            # Tool specific info: sniffles
            if tool == "sniffles":
                depth = info["DP"]
                allele_freq = info["AF"]
                allele_freq = round(float(allele_freq), 4)

            # Tool specific info: sniffles
            if tool == "snpeff":
                annotation = info["ANN"]

            # Construct entry
            variant = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "type": var_type,
                "qual": qual,
                "tool": tool,
                "depth": depth,
                "depth_fwd": depth_fwd,
                "depth_rev": depth_rev,
                "depth_ref": depth_ref,
                "depth_alt": depth_alt,
                "depth_ref_fwd": depth_ref_fwd,
                "depth_ref_rev": depth_ref_rev,
                "depth_alt_fwd": depth_alt_fwd,
                "depth_alt_rev": depth_alt_rev,
                "allele_freq": allele_freq,
                "annotation": annotation,
            }

            # Add to variants
            var_id = f"{chrom}_{pos}_{ref}_{alt}"
            if var_id not in variants:
                variants[var_id] = []
            variants[var_id].append(variant)

    # Create header
    header = [
        "chrom",
        "position",
        "type",
        "ref",
        "alt",
    ]

    # Create qual section
    for tool in tool_names:
        if tool != "snpeff":
            header += [f"{tool}_qual"]

    # Create depth section
    for tool in tool_names:
        if tool != "snpeff":
            header += [f"{tool}_depth"]

    # Add clair3 specific fields
    if "clair3" in tool_names:
        header += [
            "clair3_allele_freq",
            "clair3_depth_ref",
            "clair3_depth_alt",
        ]

    # Add lofreq specific fields
    if "lofreq" in tool_names:
        header += ["lofreq_allele_freq", "lofreq_depth_ref_fwd", "lofreq_depth_ref_rev", "lofreq_depth_alt_fwd", "lofreq_depth_alt_rev"]

    # Add freebayes specific fields
    if "freebayes" in tool_names:
        header += [
            "freebayes_allele_freq",
        ]

    # Add sniffles specific fields
    if "sniffles" in tool_names:
        header += ["sniffles_allele_freq"]

    # Add snpeff specific fields
    if "snpeff" in tool_names:
        header += [
            "annotation_allele",
            "annotation",
            "annotation_impact",
            "gene_name",
            "gene_id",
            "feature_type",
            "feature_id",
            "transcript_biotype",
            "rank",
            "hgvs_c",
            "hgvs_p",
            "cdna_pos/cdna_length",
            "cds_pos/cds_length",
            "AA.pos/AA.length",
            "distance",
            "errors_warnings_info",
        ]

    # Create blank default row entry from header
    default_row = ["" for _ in header]

    # Process variants into data table
    processed_variants = []
    for var_list in variants.values():
        # Init new row
        new_row = default_row.copy()

        # Fill consensus fields from first tool whatever it is
        curr_var = var_list[0]
        new_row[0] = curr_var["chrom"]
        new_row[1] = curr_var["pos"]
        new_row[2] = curr_var["type"]
        new_row[3] = curr_var["ref"]
        new_row[4] = curr_var["alt"]

        # Process qual section
        for tool in tool_names:
            if tool == "snpeff":
                continue

            # Check if tool exists in this var
            tool_var = [v for v in var_list if v["tool"] == tool]
            if tool_var:
                tool_var = tool_var[0]
                # Find qual position and assign
                qual_pos = header.index(f"{tool}_qual")
                new_row[qual_pos] = tool_var["qual"]

        # Process depth section
        for tool in tool_names:
            if tool == "snpeff":
                continue

            # Check if tool exists in this var
            tool_var = [v for v in var_list if v["tool"] == tool]
            if tool_var:
                tool_var = tool_var[0]
                # Find depth position and assign
                depth_pos = header.index(f"{tool}_depth")
                new_row[depth_pos] = tool_var["depth"]

        # Process clair3 specific fields
        if "clair3" in tool_names:
            clair3_var = [v for v in var_list if v["tool"] == "clair3"]
            if clair3_var:
                clair3_var = clair3_var[0]
                # Find clair3 specific positions and assign
                depth_ref_pos = header.index("clair3_depth_ref")
                depth_alt_pos = header.index("clair3_depth_alt")
                allele_freq_pos = header.index("clair3_allele_freq")
                new_row[depth_ref_pos] = clair3_var["depth_ref"]
                new_row[depth_alt_pos] = clair3_var["depth_alt"]
                new_row[allele_freq_pos] = clair3_var["allele_freq"]

        # Process lofreq specific fields
        if "lofreq" in tool_names:
            lofreq_var = [v for v in var_list if v["tool"] == "lofreq"]
            if lofreq_var:
                lofreq_var = lofreq_var[0]
                # Find lofreq specific positions and assign
                allele_freq_pos = header.index("lofreq_allele_freq")
                depth_ref_fwd_pos = header.index("lofreq_depth_ref_fwd")
                depth_ref_rev_pos = header.index("lofreq_depth_ref_rev")
                depth_alt_fwd_pos = header.index("lofreq_depth_alt_fwd")
                depth_alt_rev_pos = header.index("lofreq_depth_alt_rev")
                new_row[allele_freq_pos] = lofreq_var["allele_freq"]
                new_row[depth_ref_fwd_pos] = lofreq_var["depth_ref_fwd"]
                new_row[depth_ref_rev_pos] = lofreq_var["depth_ref_rev"]
                new_row[depth_alt_fwd_pos] = lofreq_var["depth_alt_fwd"]
                new_row[depth_alt_rev_pos] = lofreq_var["depth_alt_rev"]

        # Process freebayes specific fields
        if "freebayes" in tool_names:
            freebayes_var = [v for v in var_list if v["tool"] == "freebayes"]
            if freebayes_var:
                freebayes_var = freebayes_var[0]
                allele_freq_pos = header.index("freebayes_allele_freq")
                new_row[allele_freq_pos] = freebayes_var["allele_freq"]

        # Process sniffles specific fields
        if "sniffles" in tool_names:
            sniffles_var = [v for v in var_list if v["tool"] == "sniffles"]
            if sniffles_var:
                sniffles_var = sniffles_var[0]
                # Find sniffles specific positions and assign
                allele_freq_pos = header.index("sniffles_allele_freq")
                new_row[allele_freq_pos] = sniffles_var["allele_freq"]

        # Process snpeff specific fields
        if "snpeff" in tool_names:
            snpeff_var = [v for v in var_list if v["tool"] == "snpeff"]
            if snpeff_var:
                snpeff_var = snpeff_var[0]
                # Find snpeff specific positions and assign
                annotation_pos = header.index("annotation_allele")
                text_ann = snpeff_var["annotation"]

                # Split annotation into parts and assign to columns
                for i, allele in enumerate(text_ann):
                    parts = allele.split("|")
                    for j, part in enumerate(parts):
                        if i == 0:
                            new_row[annotation_pos + j] = part
                        else:
                            new_row[annotation_pos + j] += f" | {part}"

        # Add to processed variants
        processed_variants.append(new_row)

    # Write to file
    if output_file:
        with open(output_file, "w", encoding="UTF-8") as f:
            f.write("\t".join(header) + "\n")
            for variant in processed_variants:
                # write nothing if field is none
                f.write("\t".join([str(v) if v is not None else "" for v in variant]) + "\n")

    return variants, header, processed_variants
