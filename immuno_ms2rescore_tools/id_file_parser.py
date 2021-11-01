"""Extract PeptideRecord from identification file."""


import pandas as pd
from pyteomics import mzid
import os
import re
from typing import Dict, List
import click
from tqdm import tqdm
import tomlkit


class _IdFileParser:
    """ General file handeling methods."""

    def __init__(self, path_to_id_file, search_engine) -> None:
        self.path_to_id_file = path_to_id_file
        self.search_engine = search_engine
        self._id_df = pd.DataFrame()
        self.peprec = None

    @staticmethod
    def get_id_file_parser(search_engine, *args, **kwargs):
        parser_map = {
            "Maxquant": MaxquantFileParser,
            "Comet": CometFileParser,
            "Peaks": PeaksFileParser,
            "SpectrumMill": SpectrumMillFileParser,
            "ProteomeDiscoverer": ProteomeDiscoverer,
        }
        try:
            parser = parser_map[search_engine]
        except KeyError:
            raise NotImplementedError(
                f"Search engine '{search_engine}' not implemented."
            )
        return parser(*args, **kwargs)

    def validate_path(self):
        """ Validate the filepath."""

        if os.path.isfile(self.path_to_id_file):
            return True
        else:
            return False

    def read_to_dataframe(self, separator=None, mzid=False):
        """
        Read the id file into dataframe

        Parameters
        ---------
        mzid: bool
            True if id file is mzid, default = False
        """
        if self.validate_path():
            if mzid:
                self._id_df = mzid.DataFrame(self.path_to_id_file)
            elif (not mzid) and separator:
                self._id_df = pd.read_table(self.path_to_id_file, sep=separator)
            elif not mzid and (not separator):
                raise ValueError("Separator is required if mzid is False")

        elif self.validate_path() is False:
            raise TypeError("Not a valid filepath.")

    def peprec_to_file(self, filename: str, separator: str):
        """ Write peprec dataframe to file"""

        if not os.path.exists(filename + ".peprec"):
            self.peprec.to_csv(filename + ".peprec", separator, index=False, mode="w")
        elif os.path.exists(filename + ".peprec"):
            self.peprec.to_csv(
                filename + ".peprec", separator, index=False, header=False, mode="a"
            )

    def _load_toml(self, config):
        toml_file = ""
        with open(config, "rt") as f_in:
            for line in f_in:
                toml_file += line
        self.config = tomlkit.loads(toml_file)


class MaxquantFileParser(_IdFileParser):
    """ Parse Maxquant identification file to Peptiderecord"""

    def __init__(self, path_to_id_file, config=None) -> None:
        super().__init__(path_to_id_file, search_engine="Maxquant")
        if config:
            self._load_toml(config)
        else:
            self.config = None

    @staticmethod
    def _get_peprec_modifications(
        sequence, modification_map=None, fixed_modifications=None
    ):
        """
        Parse PEPREC modification string out of msms.txt modified sequence.

        Parameters
        ----------
        sequence: str
            msms.txt-style modified sequence
        fixed_modification: dict
            dictionary with fixed modifications
        """
        if not modification_map:
            modification_mapping = {
                "ox": "Oxidation",
                "ph": "Phospho",
                "ac": "Acetyl",
                "ca": "Carbamidomethyl",
                "cy": "Cysteinyl",
                "gl": "Gln->pyro-Glu",
                "de": "Deamidated",
            }
        else:
            modification_mapping = modification_map

        if fixed_modifications:
            inv_mod_mapping = {v: k for k, v in modification_mapping.items()}
            for mod_name, aa in fixed_modifications.items():
                mod_tag = inv_mod_mapping[mod_name]
                sequence = sequence.replace(aa, f"{aa}({mod_tag})")

        modifications = []
        if "(" not in sequence:
            modifications = "-"
        else:
            while re.match(r".*\(([^)]*)\).*", sequence):
                x = re.search(r"\(([^)]*)\)", sequence)
                modifications.append(
                    [str(x.start() - 1), modification_mapping[x.group().strip("()")]]
                )
                sequence = re.sub(r"\(([^)]*)\)", "", sequence, 1)
            modifications = ["|".join(tups) for tups in modifications]
            modifications = "|".join(modifications)
        return modifications

    def to_peprec(self):
        peprec = pd.DataFrame(
            columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "psm_score",
                "observed_retention_time",
                "Label",
                "Raw file",
            ]
        )

        if self._id_df.empty:
            self._id_df = pd.read_table(self.path_to_id_file, sep="\t")

        peprec["spec_id"] = "controllerType=0 controllerNumber=1 scan=" + self._id_df[
            "Scan number"
        ].astype(str)
        peprec["peptide"] = self._id_df["Sequence"]
        if not self.config:
            peprec["modifications"] = self._id_df["Modified sequence"].apply(
                lambda x: self._get_peprec_modifications(x)
            )
        else:
            peprec["modifications"] = self._id_df["Modified sequence"].apply(
                lambda x: self._get_peprec_modifications(
                    x,
                    self.config["modification_map"],
                    self.config["fixed_modifications"],
                )
            )
        peprec["charge"] = self._id_df["Charge"]
        peprec["psm_score"] = self._id_df["Score"]
        peprec["observed_retention_time"] = self._id_df["Retention time"]
        peprec["Label"] = self._id_df["Reverse"].isna().apply(lambda x: 1 if x else -1)
        peprec["Raw file"] = self._id_df["Raw file"]

        self.peprec = peprec


class CometFileParser(_IdFileParser):
    """ Parse Comet identification files to PeptideRecord"""

    def __init__(self, path_to_id_file) -> None:
        super().__init__(path_to_id_file, search_engine="Comet")

    @staticmethod
    def _get_peprec_modifications(sequences, mods_requiring_suffix=None):
        """Get peprec-formatted modifications."""

        if not mods_requiring_suffix:
            mods_requiring_suffix = []
        parsed_modifications = []

        for sequence in sequences:
            mod = []
            if "(" not in sequence:
                parsed_modifications.append("-")
            else:
                while re.match(r".*\(([^)]*)\).*", sequence):
                    x = re.search(r"\(([^)]*)\)", sequence)
                    loc = int(x.start())
                    name = x.group().strip("()")
                    if name in mods_requiring_suffix:
                        name = name + sequence[loc - 1]
                    mod.extend([str(loc), name])
                    sequence = re.sub(r"\(([^)]*)\)", "", sequence, 1)
                mod = "|".join(mod)
                parsed_modifications.append(mod)
        return parsed_modifications

    def to_peprec(self):
        peprec = pd.DataFrame(
            columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "psm_score",
                "observed_retention_time",
                "Label",
                "Raw file",
            ]
        )

        if self._id_df.empty:
            self._id_df = pd.read_table(self.path_to_id_file, sep="\t")

        peprec["spec_id"] = "controllerType=0 controllerNumber=1 scan=" + self._id_df[
            "ScanNr"
        ].astype(str)
        peprec["peptide"] = self._id_df["Sequence"]
        peprec["modifications"] = self._get_peprec_modifications(self._id_df["Peptide"])
        peprec["charge"] = self._id_df["Charge"]
        peprec["psm_score"] = self._id_df["Comet.SpScore"]
        peprec["observed_retention_time"] = self._id_df["RT"]
        peprec["Label"] = self._id_df["IsDecoy"].apply(lambda x: 1 if x else -1)
        raw_files = []
        for i in self._id_df.Spectrum:
            raw_file, _, _ = i.partition(".")
            raw_files.append(raw_file)
        peprec["Raw file"] = raw_files
        self.peprec = peprec


class PeaksFileParser(_IdFileParser):
    def __init__(self, path_to_id_file) -> None:
        super().__init__(path_to_id_file, search_engine="Peaks")

    @staticmethod
    def _get_peprec_modifications(modifications: List):
        """ get peprec modifications out of the peaks id file"""
        suffix_list = ["Phospho"]
        if isinstance(modifications, List):
            mods = []
            for m in modifications:
                if m["name"] in suffix_list:
                    mods.append(f"{m['location']}|{m['name'] + m['residues'][0]}")
                elif "residues" not in m.keys() and m["location"] != 0:
                    mods.append(f"-1|{m['name']}")
                else:
                    mods.append(f"{m['location']}|{m['name']}")
            return "|".join(mods)
        else:
            raise TypeError

    def convert_flatten(self, d, parent_key="", sep="_"):
        items = []
        for k, v in d.items():
            new_key = parent_key + sep + k if parent_key else k

            if isinstance(v, Dict):
                items.extend(self.convert_flatten(v, new_key, sep=sep))
            elif isinstance(v, List):
                if isinstance(v[0], Dict) and k != "Modification":
                    items.extend(self.convert_flatten(v[0], new_key, sep=sep))
                elif isinstance(v[0], Dict) and k == "Modification":
                    items.append((new_key, v))
                else:
                    items.append((new_key, v[0]))
            else:
                items.append((new_key, v))
        return items

    def to_peprec(self):

        self.peprec = pd.DataFrame(
            columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "psm_score",
                "Label",
                "Raw file",
            ]
        )
        with mzid.read(self.path_to_id_file) as reader:
            for spectrum_identification_result in tqdm(reader):
                retrieved_data = pd.Series(
                    index=[
                        "spec_id",
                        "peptide",
                        "modifications",
                        "charge",
                        "protein_list",
                        "psm_score",
                        "Label",
                        "Raw file",
                    ]
                )
                flat_dict = dict(self.convert_flatten(spectrum_identification_result))
                spec_id = (
                    flat_dict["location"]
                    .rsplit("/", 1)[1]
                    .split(".", 1)[0]
                    .replace(",", "_", 1)
                    + ":"
                    + flat_dict["spectrumID"]
                )
                retrieved_data["spec_id"] = spec_id
                retrieved_data["peptide"] = flat_dict[
                    "SpectrumIdentificationItem_PeptideSequence"
                ]
                try:
                    retrieved_data["modifications"] = self._get_peprec_modifications(
                        flat_dict["SpectrumIdentificationItem_Modification"]
                    )
                except KeyError:
                    retrieved_data["modifications"] = "-"
                retrieved_data["charge"] = flat_dict[
                    "SpectrumIdentificationItem_chargeState"
                ]
                retrieved_data["protein_list"] = flat_dict[
                    "SpectrumIdentificationItem_PeptideEvidenceRef_accession"
                ]
                retrieved_data["psm_score"] = flat_dict[
                    "SpectrumIdentificationItem_PEAKS:peptideScore"
                ]
                retrieved_data["Label"] = flat_dict[
                    "SpectrumIdentificationItem_PeptideEvidenceRef_isDecoy"
                ]
                retrieved_data["Raw file"] = (
                    flat_dict["location"]
                    .rsplit("/", 1)[1]
                    .split(".", 1)[0]
                    .replace(",", "_", 1)
                )
                self.peprec = self.peprec.append(retrieved_data, ignore_index=True)
            self.peprec["Label"] = self.peprec["Label"].apply(lambda x: -1 if x else 1)


class SpectrumMillFileParser(_IdFileParser):
    """ Parse Spectrum Mill identification file to PeptideRecord"""

    def __init__(self, path_to_id_file) -> None:
        super().__init__(path_to_id_file, search_engine="SpectrumMill")

    @staticmethod
    def _get_peprec_modifications(sequence: str):

        modification_map = {
            "m": ("M", "Oxidation"),
            "q": ("Q", "Glu->pyro-Glu"),
            "n": ("N", "Deamidated"),
            "y": ("Y", "PhosphoY"),
            "t": ("T", "PhosphoT"),
            "s": ("S", "PhosphoS"),
        }
        modification = []
        if sequence.isupper():
            return sequence, "-"
        else:
            while re.match(".*[a-z].*", sequence):
                lc = re.search("[a-z]", sequence)
                modification.append(
                    [str(lc.start() + 1), modification_map[lc.group()][1]]
                )
                sequence = re.sub("[a-z]", modification_map[lc.group()][0], sequence, 1)
            modification = ["|".join(tups) for tups in modification]
            modification = "|".join(modification)
        return sequence, modification

    @staticmethod
    def _get_filename_and_scannumber(filename: str):
        rawfilename, _, scannumber = filename.partition(".")
        scannumber = scannumber.split(".")
        scannumber = scannumber[0]
        return rawfilename, scannumber

    def to_peprec(self):
        peprec = pd.DataFrame(
            columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "protein_list",
                "psm_score",
                "observed_retention_time",
                # "Label",
                "Raw file",
            ]
        )
        if self._id_df.empty:
            self._id_df = pd.read_table(self.path_to_id_file, sep=";")

        peprec[["Raw file", "spec_id"]] = pd.DataFrame(
            self._id_df["filename"].apply(self._get_filename_and_scannumber).tolist(),
            index=self._id_df.index,
        )
        peprec[["peptide", "modifications"]] = pd.DataFrame(
            self._id_df.sequence.apply(self._get_peprec_modifications).tolist(),
            index=peprec.index,
        )
        peprec["charge"] = self._id_df["parent_charge"]
        peprec["psm_score"] = self._id_df["score"]
        peprec["observed_retention_time"] = self._id_df["retentionTimeMin"]

        self.peprec = peprec


class ProteomeDiscoverer(_IdFileParser):
    def __init__(self, path_to_id_file) -> None:
        super().__init__(path_to_id_file, search_engine="Peaks")

    @staticmethod
    def _get_peprec_modifications(modifications: List):
        """ get peprec modifications out of the proteome discoverer text, id file"""

        suffix_list = ["Phospho", "TMT6plex", "TMT12plex"]
        term_mapper = {"N-Term": 0, "C-Term": -1}
        if "(Prot)" in modifications:
            modifications = modifications.replace("(Prot)", "")
        mod_list = modifications.split("; ")
        peprec_modifications = []

        for m in mod_list:
            mod = re.search(r"\((.*?)\)", m)
            aa_loc = re.match(r"([A-Z])(\d+)", m)
            if not aa_loc:
                loc = term_mapper[m[: mod.start()]]
            else:
                aa = aa_loc.group(1)
                loc = aa_loc.group(2)

            if mod.group(1) in suffix_list:
                if loc == 0:
                    peprec_modifications.append(f"{loc}|{mod.group(1)}N")
                else:
                    peprec_modifications.append(f"{loc}|{mod.group(1)}{aa}")
            else:
                peprec_modifications.append(f"{loc}|{mod.group(1)}")
        return "|".join(peprec_modifications)

    def to_peprec(self):
        peprec = pd.DataFrame(
            columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                # "protein_list",
                "psm_score",
                "observed_retention_time",
                "q-value",
                "Raw file",
            ]
        )
        if self._id_df.empty:
            self._id_df = pd.read_table(
                self.path_to_id_file,
            )

        peprec["Raw file"] = self._id_df["Spectrum File"].apply(lambda x: x.split('.')[0])
        peprec["spec_id"] = "controllerType=0 controllerNumber=1 scan=" + self._id_df["PSMs Peptide ID"].astype(str)
        if "Sequence" in self._id_df.keys():
            peprec["peptide"] = self._id_df["Sequence"].str.upper()
        elif "Annotated Sequence" in self._id_df.keys():
            peprec["peptide"] = (
                self._id_df["Annotated Sequence"]
                .apply(lambda x: re.search(r"\].([A-z]*).\[", x).group(1))
                .str.upper()
            )
        peprec["modifications"] = self._id_df["Modifications"].apply(
            self._get_peprec_modifications
        )
        peprec["charge"] = self._id_df["Charge"]
        peprec["psm_score"] = self._id_df["DeltaScore"]
        peprec["observed_retention_time"] = self._id_df["RT [min]"]
        peprec["q-value"] = self._id_df["Percolator q-Value"]

        self.peprec = peprec


@click.command()
@click.option("--id_file", help="Filepath to identification file.")
@click.option(
    "--search_engine", help="Search engine that provided identification file."
)
@click.option("--output_filename", help="Name and directory for output peprec.")
@click.option(
    "--config",
    default=None,
    help="optional config file with fixed modifications and modification mapping",
)
def main(id_file, search_engine, output_filename, config):
    if config:
        id_file_parser = _IdFileParser.get_id_file_parser(
            search_engine, id_file, config
        )
    else:
        id_file_parser = _IdFileParser.get_id_file_parser(search_engine, id_file)
    id_file_parser.to_peprec()
    id_file_parser.peprec_to_file(output_filename, separator=" ")


if __name__ == "__main__":
    main()
