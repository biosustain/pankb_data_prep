import pandas as pd
import json
import argparse
from bs4 import BeautifulSoup
from bs4.element import Tag, NavigableString, Comment
from urllib.request import urlopen
from urllib.parse import urlparse, parse_qs

IMODULONDB_URL = "https://imodulondb.org/index.html"

def initialize_parser(parser):
    parser.description = "Parse iModulonDB for general info."
    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output file.",
    )


def im_init(
    output_path
):
    output = urlopen(IMODULONDB_URL).read()
    im_html = output.decode('utf-8')

    soup = BeautifulSoup(im_html, 'html.parser')
    org_acc = soup.find(id="organismAccordion")

    org_d = {}

    last_org_name = None
    for c in org_acc.contents:
        if isinstance(c, Comment):
            last_org_name = str(c.string).strip()
        elif isinstance(c, Tag):
            if last_org_name is None:
                raise Exception()
            org_name = last_org_name
            org_short = None
            org_datasets = []
            last_org_name = None
            dl = c.find("ul", class_="nav")
            for dataset in dl.find_all("li"):
                href = dataset.a['href']
                href_parts = parse_qs(urlparse(href).query)
                if org_short is None:
                    org_short = href_parts["organism"][0]
                org_datasets.append((href_parts["dataset"][0], href))
            org_d[org_name] = {"short": org_short, "datasets": org_datasets}
    with open(output_path, "w") as f:
        json.dump(org_d, f)

def run(args):
    im_init(
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
