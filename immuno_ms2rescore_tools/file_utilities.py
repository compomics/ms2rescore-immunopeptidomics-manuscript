"""utilities for handeling different filetypes"""

from os.path import isdir
from re import escape
from numpy.core.numeric import full
import pandas as pd
import os
from tqdm import tqdm
from pyteomics.auxiliary import target_decoy
import re
from pyteomics import mgf


class FileHandeling:
    """Standard filehandeling methods"""

    def __init__(self) -> None:
        self.filelist = []
        pass

    def retrieve_files(self, path, file_extension=None, recursive=False):
        if os.path.isfile(path):
            self.filelist.append(path)
        elif os.path.isdir(path):
            dirs = []
            files = []
            for object in os.listdir(path):
                object_path = os.path.join(path, object)
                if os.path.isfile(object_path) and object_path.endswith(
                    "." + file_extension
                ):
                    files.append(object_path)
                elif os.path.isdir(object_path):
                    dirs.append(object_path)
            if recursive:
                for directories in dirs:
                    dirs = []
                    FileHandeling.retrieve_files(self, directories, file_extension)
                os.chdir("..")
            self.filelist.extend(files)
        else:
            raise TypeError("Not a valid file or directory")

    def combine_peprec_files(self, filelist=None):
        if filelist:
            self.filelist = filelist
        peprecs = []
        for filename in self.filelist:
            peprecs.append(pd.read_table(filename, sep=" "))
        peprec = pd.concat(peprecs, ignore_index=True, join="inner")
        return peprec


class MascotGenericFormat(FileHandeling):
    """Methods for MGF files"""

    def __init__(self, mgf_folder) -> None:
        super().__init__()
        self.retrieve_files(mgf_folder, file_extension="mgf")

    def count_spectra(self, verbose=True):
        if not self.filelist:
            print("No files listed")

        total = 0
        for mgf_file_path in tqdm(self.filelist):
            count = 0
            with open(mgf_file_path, "r") as f:
                for line in f:
                    if line.rstrip("\n") == "BEGIN IONS":
                        count += 1
            filename = mgf_file_path.rsplit("/", 1)[1]
            total += count
            if verbose:
                print(f"{filename}: {count}")
        print(f"Total spectra: {total}")

    def get_spec_id(self):
        spec_list = []
        for mgf_file in self.filelist:
            with open(mgf_file, "r") as f:
                for line in f:
                    if "TITLE" in line:
                        spec_list.append(line[6:].strip())

    def scan_mgf(self, peprec_in, spec_id_name, outname="scan_mgf_result.mgf", usi=False):

        with open(outname, mode="w") as out:
            for mgf_file in tqdm(self.filelist):
                if not usi:
                    spec_dict = dict(
                        (
                            re.search(r"scan(?:\=|\:)(\d+)", v).group(1),
                            peprec_in[spec_id_name][k],
                        )
                        for k, v in peprec_in[
                            (
                                peprec_in["Raw file"]
                                == mgf_file.rsplit("/", 1)[1].split(".", 1)[0]
                            )
                        ]["spec_id"]
                        .to_dict()
                        .items()
                    )
                if usi:
                    spec_dict = dict(
                        (
                            re.search(r"mzspec:(unpublished|PXD[0-9]{6}):(\S*):scan:\d+", v).group(0),
                            peprec_in[spec_id_name][k],
                        )
                        for k, v in peprec_in[
                            (
                                peprec_in["Raw file"]
                                == mgf_file.rsplit("/", 1)[1].split(".", 1)[0]
                            )
                        ]["spec_id"]
                        .to_dict()
                        .items()
                    )
                found = False
                with open(mgf_file, "r") as f:
                    for line in f:
                        if "TITLE" in line:
                            if not usi:
                                match = re.search(r"scan(?:\=|\:)(\d+)", line[6:].strip())
                                scan_id = match.group(1)
                            elif usi:
                                match = re.match(r"mzspec:(unpublished|PXD[0-9]{6}):(\S*):scan:\d+", line[6:].strip())
                                scan_id = match.group(0)

                            if scan_id in spec_dict:
                                found = True
                                out.write("BEGIN IONS\n")
                                out.write("TITLE=" + spec_dict[scan_id] + "\n")
                                continue
                        elif "END IONS" in line:
                            if found:
                                out.write(line)
                                found = False
                        elif found and line[-4:] != "0.0\n":
                            out.write(line)

    def retrieve_masses(self, psm_id):
        spectrum_dict = dict()
        found = False
        for mgf_file in self.filelist:
            for spectrum in mgf.read(mgf_file):
                if psm_id == spectrum["params"]["title"]:
                    spectrum_dict["identifier"] = spectrum["params"]["title"]
                    spectrum_dict["precursor_mz"] = spectrum["params"]["pepmass"][0]
                    spectrum_dict["precursor_charge"] = spectrum["params"]["charge"][0]
                    spectrum_dict["mz"] = spectrum["m/z array"]
                    spectrum_dict["intensity"] = spectrum["intensity array"]
                    spectrum_dict["retention_time"] = float(
                        spectrum["params"]["rtinseconds"]
                    )
                    found = True
                    break
            if found:
                break
        return spectrum_dict

    def check_mgf_file_presence(self, raw_file_list: list):
        "Check if all mgf files are present in your folder, given a raw file list"

        missing_mgf_files = []

        mgf_files = [x.rsplit("/", 1)[1].split(".", 1)[0] for x in self.filelist]
        for raw_file in raw_file_list:
            if raw_file in mgf_files:
                pass
            else:
                missing_mgf_files.append(raw_file)
        return missing_mgf_files


class PeptideRecord(FileHandeling):
    """Methods for peprec files"""

    def __init__(self, peprec) -> None:
        super().__init__()
        if isinstance(peprec, str):
            self.peprec = pd.read_table(peprec, sep=" ")
            if "/" in peprec:
                self.peprec_name = peprec.rsplit("/", 1)[1].split(".", 1)[0]
            else:
                self.peprec_name = peprec.split(".", 1)[0]
        elif isinstance(peprec, pd.DataFrame):
            self.peprec = peprec
            self.peprec_name = "peprec"

    def select_unique_peptide(self):
        self.peprec.sort_values("psm_score", ascending=False)
        self.peprec.drop_duplicates(
            ["peptide", "modifications", "charge"], keep="first", inplace=True
        )

    def select_duplicates(self):
        duplicates = self.peprec.duplicated(["peptide", "modifications", "charge"], keep=False)
        self.peprec = self.peprec[duplicates]

    def _get_double_modifications(self):
        unique_aa_modifications = dict()
        self._suffix_list = []
        for tuple in self.peprec.itertuples():
            s = re.findall(r"[\d]*.\|[A-Za-z0-9]*", tuple.modifications)
            for match in s:
                loc = int(re.search(r"(\d+)", match).group(0))
                modification = re.search(r"[A-Za-z0-9]*$", match).group()
                if loc == 0:
                    aa = "N-term"
                elif loc == len(tuple.peptide) + 1:
                    aa = "C-term"
                else:
                    aa = tuple.peptide[loc - 1]

                if (
                    modification in unique_aa_modifications.keys()
                    and aa == unique_aa_modifications[modification]
                ):
                    pass
                elif (
                    modification in unique_aa_modifications.keys()
                    and aa != unique_aa_modifications[modification]
                ):
                    if modification in self._suffix_list:
                        pass
                    elif modification not in self._suffix_list:
                        self._suffix_list.append(modification)
                else:
                    unique_aa_modifications[modification] = aa

    def add_modification_suffix(self):
        self._get_double_modifications()

        def sequence_location_fix(location: str):
            location = int(location)
            if location > 0:
                return location - 1
            else:
                return location

        if not self._suffix_list:
            pass
        else:
            print(self._suffix_list)
            modified_peprec = self.peprec[
                (self.peprec["modifications"] != "-")
                & (self.peprec["modifications"] != "")
                & (~self.peprec["modifications"].isna())
            ]
            peptide_modifications = []
            for tuple in modified_peprec.itertuples(name="peprec", index=False):
                new_mods = []
                sequence = tuple.peptide
                loc = tuple.modifications.split("|")[::2]
                mod = tuple.modifications.split("|")[1::2]
                for loc, mod in zip(loc, mod):
                    if mod in self._suffix_list:
                        new_mods.append(
                            f"{loc}|{mod + sequence[sequence_location_fix(loc)]}"
                        )
                    else:
                        new_mods.append(f"{loc}|{mod}")

                peptide_modifications.append("|".join(new_mods))
            modified_peprec = modified_peprec.drop("modifications", axis=1)

            modified_peprec.insert(2, "modifications", peptide_modifications)
            peprec = pd.concat(
                [
                    modified_peprec,
                    self.peprec[
                        ~(self.peprec["modifications"] != "-")
                        & (self.peprec["modifications"] != "")
                        & (~self.peprec["modifications"].isna())
                    ],
                ]
            ).sort_index()

            self.peprec = peprec

    def calculate_qvalues(self):
        runs = self.peprec["Raw file"].unique()
        qvalues_list = []
        for run in runs:
            subset = pd.DataFrame(
                target_decoy.qvalues(
                    self.peprec.loc[self.peprec["Raw file"] == run],
                    is_decoy=self.peprec["Label"] == -1,
                    key="psm_score",
                    reverse=True,
                    formula=1,
                    full_output=True,
                )
            )
            qvalues_list.append(subset)
        self.qvalues = pd.concat(qvalues_list)

    def filter_peprec_on_qvalue(self, limit=0.01):
        self.calculate_qvalues()
        self.qvalues = self.qvalues.loc[self.qvalues["q"] < limit]
        self.peprec = self.peprec[self.peprec.isin(self.qvalues)].dropna()

    def filter_decoys(self):
        self.peprec = self.peprec.loc[self.peprec.Label == 1]

    def count_decoys(self):
        decoys = sum(self.peprec.Label == 1)
        return decoys

    def remove_peptides_without_spectrum(self, raw_files_to_remove: list):
        if not raw_files_to_remove:
            pass
        else:
            for raw_file in raw_files_to_remove:
                self.peprec.drop(
                    self.peprec.loc[self.peprec["Raw file"] == raw_file].index,
                    inplace=True,
                )

    def create_usi(self, identifier):
        self.peprec.rename(columns={"Raw file": "Raw_file"}, inplace=True)
        usi_list = []
        if all(
            [
                True
                if re.match(r"mzspec:(unpublished|PXD[0-9]{6}):(\S*):scan:\d+", x)
                else False
                for x in self.peprec.spec_id
            ]
        ):
            usi_list = self.peprec["spec_id"]
        elif identifier:
            for usi in self.peprec.itertuples():
                usi_list.append(
                    f"mzspec:{identifier}:{usi.Raw_file}:scan:{usi.spec_id.rsplit('=', 1)[1]}"
                )
        self.peprec.rename(columns={"Raw_file": "Raw file"}, inplace=True)
        self.peprec.insert(0, "USI", usi_list)

    def peprec_to_file(self, filename: str, separator: str):
        """ Write peprec dataframe to file"""

        if not os.path.exists(filename + ".peprec"):
            self.peprec.to_csv(filename + ".peprec", separator, index=False, mode="w")
        elif os.path.exists(filename + ".peprec"):
            self.peprec.to_csv(
                filename + ".peprec", separator, index=False, header=False, mode="a"
            )

    def to_prosit_csv(self, ce):
        """Convert peprec to prosit csv for predicting spectral intensities with prosit"""

        prosit_csv = pd.DataFrame()
        prosit_csv["precursor_charge"] = self.peprec["charge"]
        prosit_csv["collision_energy"] = ce
        prosit_csv["modified_sequence"] = self.acquire_modified_seq()

        return prosit_csv

    @staticmethod
    def _get_modified_sequence(sequence_mod_tuple):
        """
        Combine peptide sequence and peprec modification into modified sequence

        input:
        sequence_mod_tuple: tuple of seq and peprec style modification

        output:
        modified sequence
        """
        inverse_modification_map = {
            "Oxidation" : "ox",
            "Phospho" : "ph",
            "Acetyl" : "ac",
            "Carbamidomethyl" : "ca",
            "Cysteinyl" : "cy",
        }
        seq, peprec_mod = sequence_mod_tuple
        if peprec_mod == "-":
            return seq

        locations = peprec_mod.split("|")[::2]
        modifications = peprec_mod.split("|")[1::2]

        for i,(loc, mod) in enumerate(zip(locations, modifications)):
            if loc == "-1":
                seq += f"({inverse_modification_map[mod]})"

            seq = seq[:int(loc)+ i*4 ] + f"({inverse_modification_map[mod]})" + seq[int(loc) + i*4:]

        return seq

    def acquire_modified_seq(self):
        """Combine peptide and modifications column of peprec to modified sequences"""

        loc_mod_tuples = pd.Series(zip(self.peprec["peptide"], self.peprec["modifications"]))
        modified_sequences = loc_mod_tuples.apply(self._get_modified_sequence)

        return modified_sequences
