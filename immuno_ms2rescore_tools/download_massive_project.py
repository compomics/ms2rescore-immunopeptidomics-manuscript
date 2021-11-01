"""
## Download massIVE Project
Download PRIDE project files for a given PRIDE identifier. With the `-f`
argument certain file types can be chosen for download.

*Download_PRIDE_Project.py*
**Input:** massIVE Archive identifier
**Output:** Downloaded files

"""

import os
import argparse
import socket


from ftplib import FTP
import ftplib
from time import sleep
from tqdm import tqdm
import pandas as pd


""" Global variables """

MY_DIRS = []  # global
MY_FILES = []  # global
CURDIR = ""  # global


def argument_parser():
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description="Download files from massIVE database for a given identifier "
    )
    parser.add_argument(
        "massive_identifier",
        action="store",
        help="Specific identiefier for the MassIVE database",
    )
    parser.add_argument(
        "-f",
        dest="filetypes",
        action="store",
        help="filetypes to download from massIVE library",
        default="raw",
    )
    parser.add_argument(
        "-s",
        dest="store",
        action="store",
        default=os.getcwd(),
        help="location to store the files, if not givin files will be stored in current directory",
    )
    parser.add_argument(
        "-n",
        action="store_true",
        dest="directory",
        help=" If argument is given new directory will be made to store the file with name of massIVE identifier",
    )
    """parser.add_argument(
        "-t",
        dest= "timeout",
        action= 'store',
        help= "The amount of time before the program stops",
        default= "2000",
        type= int,
    )"""
    args = parser.parse_args()
    return args


def get_dirs(ln):
    """ list directories and files with certain filetypes"""

    args = argument_parser()
    global MY_DIRS
    global MY_FILES
    cols = ln.split(" ")
    objname = cols[len(cols) - 1]  # file or directory name
    if ln.startswith("d"):
        MY_DIRS.append(objname)
    else:
        if objname.endswith("." + args.filetypes):
            MY_FILES.append(os.path.join(CURDIR, objname))  # full path
    return


def check_dir(adir, ftp):
    """
    Search all subdirectories on a ftp connection given a starting directory

    Parameters
    ----------
    adir: str
        directory to start searching
    ftp: object
        established ftp connection
    """
    global MY_DIRS
    global MY_FILES  # let it accrue, then fetch them all later
    global CURDIR

    MY_DIRS = []
    gotdirs = []  # local
    CURDIR = ftp.pwd()

    ftp.cwd(adir)
    CURDIR = ftp.pwd()
    ftp.retrlines("LIST", get_dirs)
    gotdirs = MY_DIRS
    sleep(1)
    for subdir in gotdirs:
        MY_DIRS = []
        check_dir(subdir, ftp)  # recurse
    ftp.cwd("..")  # back up a directory when done here
    return


def download_file(directory, file_name):
    """
    Download the file given a directory

    Parameters
    ---------
    directory: str
        file directory on ftp site
    file_name: str
        name for the local file
    """

    with FTP("massive.ucsd.edu", user="anonymous", passwd="anonymous@") as ftp:
        ftp.sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
        ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, 75)
        ftp.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, 60)
        with open(file_name, "wb") as rf:
            total = ftp.size(directory)
            with tqdm(total=total, unit_scale=True, unit_divisor=1024) as pbar:

                def cb(data):
                    pbar.update(len(data))
                    rf.write(data)

                ftp.retrbinary("RETR {}".format(directory), cb)


def main():
    args = argument_parser()
    with FTP("massive.ucsd.edu", user="anonymous", passwd="anonymous@") as ftp:

        print("getting files...")
        check_dir(args.massive_identifier, ftp)  # directory to start in
        print("Total files found: {} ".format(len(MY_FILES)))
        ftp.cwd("/.")  # change to root directory for downloading

        os.chdir(args.store)

        if args.directory:
            if not os.path.exists(os.path.join(args.store, args.massive_identifier)):
                os.mkdir(
                    os.path.join(args.store, args.massive_identifier)
                )  # make a directory to save the files
            os.chdir(os.path.join(args.store, args.massive_identifier))

        existing_files = os.listdir()
        failed_files = []

    for f in tqdm(MY_FILES):
        dir, sep, filename = f.rpartition("/")

        if filename in existing_files:
            continue
        else:
            try:
                download_file(f, filename)
            except ftplib.error_temp:
                failed_files.append(filename)
                continue

    print("Following files failed to download:")
    print(failed_files)


if __name__ == "__main__":
    main()
