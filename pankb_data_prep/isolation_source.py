import pandas as pd
import requests
from .utilities import get_genome_list
import argparse


def initialize_parser(parser):
    parser.description = "Fetch the source of the isolates from NCBI."
    parser.add_argument(
        "--genomes",
        "-g",
        type=str,
        required=True,
        help="Genome IDs to process.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )


def get_biosample_info(assembly_id):
    """
    Fetch biosample info from NCBI Datasets API v2alpha.
    Returns a dict with biosample_accession, isolation_source, geo_loc_name, and Country.
    """
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_id}/dataset_report"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code != 200:
            return {
                "biosample_accession": "Missing",
                "isolation_source": "Failed to retrieve the data",
                "geo_loc_name": "Failed to retrieve the data",
                "Country": "Failed to retrieve the data",
            }

        data = response.json()
        reports = data.get("reports", [])
        if not reports:
            return {
                "biosample_accession": "Missing",
                "isolation_source": "Missing",
                "geo_loc_name": "Missing",
                "Country": "Missing",
            }

        biosample = reports[0].get("assembly_info", {}).get("biosample", {})
        biosample_accession = biosample.get("accession", "Missing")
        isolation_source = biosample.get("isolation_source", "Missing")
        geo_loc_name = biosample.get("geo_loc_name", "Missing")

        # Extract country from geo_loc_name (format: "Country:Region:City")
        if geo_loc_name and geo_loc_name != "Missing":
            country = geo_loc_name.split(":")[0]
        else:
            country = "Missing"

        return {
            "biosample_accession": biosample_accession,
            "isolation_source": isolation_source,
            "geo_loc_name": geo_loc_name,
            "Country": country,
        }

    except Exception as e:
        return {
            "biosample_accession": "Missing",
            "isolation_source": "Failed to retrieve the data",
            "geo_loc_name": "Failed to retrieve the data",
            "Country": "Failed to retrieve the data",
        }


def find_isolation_source(genome_ids, isolation_source_path):
    results = []
    for genome_id in genome_ids:
        info = get_biosample_info(genome_id)
        info["genome_id"] = genome_id
        results.append(info)

    samples = pd.DataFrame(results)
    # Reorder columns to match original output format
    samples = samples[["genome_id", "biosample_accession", "isolation_source", "geo_loc_name", "Country"]]
    samples.to_csv(isolation_source_path, index=False)


def run(args):
    genomes = get_genome_list(args)
    find_isolation_source(genomes, args.output)


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
