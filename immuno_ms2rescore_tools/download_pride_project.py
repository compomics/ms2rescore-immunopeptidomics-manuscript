"""
## Download PRIDE Project
Download PRIDE project files for a given PRIDE identifier. With the `-f`
argument certain file types can be chosen for download.

*Download_PRIDE_Project.py*
**Input:** PRIDE Archive identifier
**Output:** Downloaded files, sorted in folders by file type
"""

import time
import json
import os
import requests
import argparse
import re
import socket
import ftplib
from collections import defaultdict
from tqdm import tqdm

import pandas as pd


def argument_parser():
    parser = argparse.ArgumentParser(
        description="Download files from PRIDE Archive for a given project."
    )
    parser.add_argument(
        "pxd_identifier",
        action="store",
        help="PXD identifier of project from which to download files",
    )
    parser.add_argument(
        "-p",
        dest="pattern",
        action="store",
        help="Regex pattern matching to files to be downloaded",
    )
    parser.add_argument(
        "-f",
        dest="filetype",
        action="store",
        help="filetype to download (msf, raw, txt, zip...)",
    )
    parser.add_argument(
        "-m",
        action="store_true",
        dest="metadata",
        help="optional paramater to download metadata",
    )
    args = parser.parse_args()
    return args


def check_pxd_id(pxd_identifier):
    """
    Assert if the project data for a given PXD identifier is accessable through
    the PRIDE Archive API.
    """
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{pxd_identifier}"
    response = requests.get(url)
    status_codes = {
        "401": "Unauthorized",
        "403": "Forbidden",
        "404": "Not Found",
        "500": "Internal server error",
        "400": "Bad request",
    }
    assert (
        response.status_code == 200
    ), f"Error code{str(response.status_code)} :\
         {status_codes[str(response.status_code)]}"


def get_files_df(pxd_identifier, filetype=None, pattern=None):

    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={pxd_identifier}"
    check_pxd_id(pxd_identifier)

    response = requests.get(url).json()
    files_list = defaultdict(list)

    for files in response:
        for file_locations in files["publicFileLocations"]:
            if file_locations["name"] == "FTP Protocol":
                files_list["filename"].append(files["fileName"])
                files_list["ftp"].append(file_locations["value"])
            else:
                pass
    filelist = pd.DataFrame(dict(files_list))
    filelist["extension"] = filelist["filename"].apply(lambda x: x.rsplit(".", 1)[1])

    if pattern:
        filelist = filelist[filelist["filename"].str.contains(pattern)]
    if filetype:
        filelist = filelist[filelist["extension"] == filetype]

    return filelist


def FTP_server_download(ftp_link):
    filename = ftp_link.rsplit("/", 1)[1]
    directory = re.search(r"\/pride[\S]*", ftp_link).group(0)

    with ftplib.FTP("ftp.pride.ebi.ac.uk", user="anonymous", passwd="anonymous@") as ftp:
        ftp.sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
        ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, 75)
        ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, 60)
        with open(filename, "wb") as rf:
            ftp.retrbinary("RETR {}".format(directory), rf.write)


def check_present_files(path):
    present_files = []
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            present_files.append(file)

    return present_files


def run():
    args = argument_parser()

    if args.metadata:
        print("Downloading meta data...")
        url = (
            f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{args.pxd_identifier}"
        )
        response = requests.get(url).json()
        with open("pxd_project_metadata.json", "w") as f:
            f.write(json.dumps(response, indent=4))

    # Download files
    print("Downloading files...")
    file_df = get_files_df(args.pxd_identifier, args.filetype, args.pattern)
    filelist = list(file_df["ftp"])
    files_in_dir = check_present_files(os.getcwd())
    retry_files = []
    count = 0
    for file in tqdm(filelist):
        if re.search(r"[^\/]*." + args.filetype, file).group(0) in files_in_dir:
            continue
        count += 1
        try:
            FTP_server_download(file)
        except ftplib.error_perm:
            retry_files.append(file)
            continue
        if count == 30:
            time.sleep(60)
            count = 0

    failed_files = []
    count = 0
    for file in retry_files:
        count += 1
        try:
            FTP_server_download(file)
        except ftplib.error_perm:
            retry_files.append(file)
            continue
        if count == 30:
            time.sleep(60)
            count = 0
    if failed_files:
        print(f"Following files failed to download:{failed_files}")


if __name__ == "__main__":
    run()
