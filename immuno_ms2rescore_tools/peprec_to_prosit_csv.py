import pandas as pd
import click

from immuno_ms2rescore_tools import file_utilities

@click.command()
@click.option("--peprec", help="peprec input file")
@click.option("--ce", help="collision energy used")
@click.option("--output_name", help="output_filename")
def main(peprec, ce, output_filename):
    peprec = file_utilities.PeptideRecord(peprec)
    prosit_csv = peprec.to_prosit_csv(ce)
    prosit_csv.to_csv(f"{output_filename}.csv", sep=",")

if __name__ == "__main__":
    main()

