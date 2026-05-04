"""Module for running twirl to re-solve WCS headers 
in a speciifc LCO-Sinistro photometry .fits files. 
(files must first be downloaded from the 
LCO Archive and unzipped using `funpack`).

NOTE: This module requires an internet connection to access Gaia DR3.
And also requires that it be called via command line with the path to the
directory only. The directory should only contain .fits files to be re-solved."""

import os
import sys
import warnings
from datetime import datetime
from glob import glob
from requests import HTTPError
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import pytz
from twirl import compute_wcs, gaia_radecs, find_peaks
from twirl.geometry import sparsify
from .resourcefiles import load_wcs_keys

wcs_keys = load_wcs_keys()

def resolve_wcs(file, backup=False):
    """ Re-solves an LCO Sinistro WCS in the file header. 
    Takes in a .fits file and returns the same file, with a newly re-solved WCS using twirl's 
    interface for finding sources and matching to Gaia DR3.
    Args:
        file (.fits): .fits file from LCO Archive
        backup (bool, optional): Whether to create a backup of the original file.
    """
    assert os.path.isfile(file), f"{file} is not a valid file"
    assert file.endswith('.fits'), f"{file} is not a .fits file"

    with fits.open(file, mode='update', save_backup=backup) as hdu:
        if (hdu[0].shape != (4096, 4096)):
            hdunum = 1
        else:
            hdunum = 0
        assert hdu[hdunum].header['EXPTIME']>0, f"{file} has zero exposure time" # pylint: disable=E1101

        hdu1 = hdu[hdunum]
        if hdu1.header['EXPTIME']<1: # pylint: disable=E1101
            print(f"Exposure <1 seconds, so file {file} should be skipped", file=sys.stdout)
            os.rename(file, os.path.join(os.path.split(file)[0],
                                        "twirl_failed", os.path.split(file)[1]))
            return
        data, original_wcs = hdu1.data, WCS(hdu1.header) # pylint: disable=E1101
        xy = find_peaks(data)[0:20]
        time = datetime.strptime(hdu1.header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f') # pylint: disable=E1101
        time_aware = pytz.utc.localize(time)
        fov = (data.shape * proj_plane_pixel_scales(original_wcs))[0]
        center = original_wcs.pixel_to_world(*np.array(data.shape) / 2)
        radecs = gaia_radecs(center, 1.2*fov, limit=50, dateobs=time_aware)
        radecs = sparsify(radecs, 0.01)
        wcs = compute_wcs(xy, radecs[0:20], tolerance=5)
        new_wcs = wcs.to_fits()
        for key in wcs_keys:
            if key=="PC1_1":
                rep = "CD1_1"
            elif key=="PC1_2":
                rep = "CD1_2"
            elif key=="PC2_1":
                rep = "CD2_1"
            elif key=="PC2_2":
                rep = "CD2_2"
            else:
                rep = str(key)
            hdu1.header[rep] = new_wcs[0].header[key] # pylint: disable=E1101
        hdu.flush()
        return
    
def run_alignment(directory, backup=False, test=False):
    """ Runs resolve_wcs on all LCO Sinistro .fits images in a given directory. 
    Also moves any files that twirl fails to re-solve the WCS for to a 
    subdirectory called "twirl_failed" at the same path as the original directory.

    NOTE: This function does **not** work with the zipped .fits files
    that can be downloaded off the LCO Science Archive. The files must all be
    unpacked beforehand (e.g. "funpack *.fz")

    Args:
        directory (str): Path to directory containing .fits files to be re-solved.
        backup (bool, optional): Whether to create a backup of the original files.
        test (bool, optional): Whether to run in test mode.
    """
    warnings.filterwarnings('ignore')
    filelist = sorted(glob(os.path.join(directory, "*.fits")))
    os.makedirs(os.path.join(directory, "twirl_failed"), exist_ok = True)
    print(filelist)
    for f in filelist:
        try:
            resolve_wcs(f, backup=backup)
        except (ValueError, HTTPError) as e:
            print(f"Error {e} so file {f} should be skipped", file=sys.stdout)
            if test is not True:
                os.rename(f, os.path.join(directory, "twirl_failed", os.path.split(f)[1]))
    if test is True:
        return filelist
