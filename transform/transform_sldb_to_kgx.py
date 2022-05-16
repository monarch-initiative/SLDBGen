#!/usr/bin/env python
# -*- coding: utf-8 -*-

from koza.cli_runner import transform_source #type: ignore


"""
Script for transforming SLDB output tsv to KGX tsv,
compatible with Biolink Model.
"""

def parse_input(data_file: str = "../SL_data.tsv") -> None:
    """
    Load SLDB data tsv.
    """
    print(f"Transforming {data_file}")

    transform_source(source="sl_data_config.yaml", 
                    output_dir=".",
                    output_format="tsv")

if __name__ == "__main__":
    parse_input()