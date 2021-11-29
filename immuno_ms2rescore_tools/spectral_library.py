"create spectral library "

import pandas as pd
from immuno_ms2rescore_tools import file_utilities
import click


class Spectrallibrary:
    def __init__(self, peprec, mgf_file_list) -> pd.DataFrame:
        self.df = file_utilities.PeptideRecord(peprec)
        self.mgf_folder = mgf_file_list

    def create_spectral_library_from_pep(self, identifier):
        "Create a spectral library peprec with concomitant mgf file"
        if "decoys" in self.df.peprec.keys():
            number_of_decoys = self.df.count_decoys()
            if number_of_decoys > 0:
                print("Filtering out q values lower than 0.01")
                self.df.filter_peprec_on_qvalue()
                print("Filtering out decoys")
                self.df.filter_decoys()
        elif "q-values" in self.df.peprec.keys():
            self.df.loc[self.df["q-values"] < 0.01]
            self.df.drop("q-values", axis=1, inplace=True)
        else:
            print("no decoys/q-values present")
        print("Selecting unique peptides")
        self.df.select_unique_peptide()
        print(f"{len(self.df.peprec.peptide)} unique peptides found")
        print("Checking if modifications are unique")
        self.df.add_modification_suffix()

        print("Gathering mgf files in folder")
        self.mgf = file_utilities.MascotGenericFormat(self.mgf_folder)
        print("Removeving peptides without spectra (mgf file not present)")
        missing_mgf = self.mgf.check_mgf_file_presence(
            self.df.peprec["Raw file"].unique()
        )
        self.df.remove_peptides_without_spectrum(missing_mgf)
        print(f"Final number unique peptides: {len(self.df.peprec)}")
        print("Gathering peptides spectra in one mgf file")
        if identifier:
            self.df.create_usi(identifier)
            self.mgf.scan_mgf(
                self.df.peprec, spec_id_name='USI', outname="spec_lib_" + self.df.peprec_name + ".mgf"
            )
            self.df.peprec = self.df.peprec.drop("spec_id", axis=1)
            self.df.peprec = self.df.peprec.rename(columns={"USI": "spec_id"})
            self.df.peprec["Raw file"] = "spec_lib_" + self.df.peprec_name
            self.df.peprec.to_csv(
                "spec_lib_" + self.df.peprec_name + ".peprec",
                sep=" ",
                index=False,
                header=True,
                mode="w",
            )
        else:
            self.mgf.scan_mgf(
                self.df.peprec, spec_id_name="spec_id", outname="spec_lib_" + self.df.peprec_name + ".mgf"
            )


@click.command()
@click.option("--peprec", help="peprec to create spectral library from")
@click.option("--mgf_folder", help="mgf folder/file with concomitant spectra")
@click.option("--identifier", default=None, help="Idenitifier for Universal Spectrum Identifier")
def main(peprec, mgf_folder, identifier):
    spectral_lib = Spectrallibrary(peprec, mgf_folder)
    spectral_lib.create_spectral_library_from_pep(identifier)


if __name__ == "__main__":
    main()
