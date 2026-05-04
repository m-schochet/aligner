import os
import pytest
import requests
from astropy.io import fits
from tests import TESTDATA
from aligner.twirler import run_alignment

def get_json(url, token=None):
    """Turn a URL into JSON for requests"""
    headers = {}
    if token and token != "AUTHTOKEN":
        headers["Authorization"] = f"Token {token}"

    r = requests.get(url, headers=headers)
    r.raise_for_status()
    return r.json()

def download_file(url, filename, chunk_size=8192):
    """Download a file from a URl turned into a JSON"""
    headers = {}
    mode = "wb"

    with requests.get(url, headers=headers, stream=True) as r:
        if r.status_code not in (200, 206):
            r.raise_for_status()

        with open(os.path.join(TESTDATA, filename), mode) as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)


def download_frame(frame_id, authtoken, frame_url):
    """Wrapper for downloading files"""
    json_url = f"{frame_url}{frame_id}/"
    metadata = get_json(json_url, authtoken)

    filename = metadata["filename"]
    url = metadata["url"]

    print(f"Downloading {filename}...")
    download_file(url, filename)

    return filename

def unpack_fz(filename):
    """Take a downloaded .fz file and unpack it to a .fits file"""
    with fits.open(filename) as hdul:
        hdul.writeto(filename.replace('.fz', ''), overwrite=True)

@pytest.fixture(scope="session")
def downloaded_frame():
    """Actual download script for a specific 1m-telescope image from the LCO Archive"""
    frames = [83071988]
    authtoken = "fa6f5b585d7914fd76c2e0047335329efc042d14"
    frame_url = "https://archive-api.lco.global/frames/"

    downloaded_files = []

    print(f"Downloading {len(frames)} frames")

    for f in frames:
        filename = download_frame(f, authtoken, frame_url)
        downloaded_files.append(filename)

    unpack_fz(os.path.join(TESTDATA, downloaded_files[0]))

    # Give files to tests
    yield downloaded_files[0].replace('.fz', '')
    test_1m_alignment(downloaded_files[0].replace('.fz', ''))

    # Cleanup
    for f in downloaded_files:
        if os.path.exists(os.path.join(TESTDATA, f)):
            print(f"Removing {f}")
            os.remove(os.path.join(TESTDATA, f))
            print(f"Removing {f.replace('.fz', '')}")
            os.remove(os.path.join(TESTDATA, f).replace('.fz', ''))

def test_1m_alignment(downloaded_frame):
    """Only test in this file"""
    file_path = downloaded_frame
    result = run_alignment(TESTDATA, backup=False, test=True)
    test_filepath = os.path.join(TESTDATA, file_path)
    with fits.open(test_filepath) as hdu:
        if (hdu[0].shape != (4096, 4096)):
            hdunum = 1
        else:
            hdunum = 0
        assert hdu[hdunum].header['EXPTIME'] > 0, f"{test_filepath} has zero exposure time"
        assert hdu[hdunum].shape == (4096, 4096), f"{test_filepath} has incorrect shape"
        assert hdu[hdunum].header['ORIGIN'] == 'LCOGT', f"{test_filepath} is not an LCO image"
        assert hdu[hdunum].header['CD1_1'] == pytest.approx(0.00010827919893101, abs=1e-8), f"{test_filepath} has incorrect CD1_1"
        assert hdu[hdunum].header['CD1_2'] == pytest.approx(3.7373047427083e-07, abs=1e-8), f"{test_filepath} has incorrect CD1_2"
        assert hdu[hdunum].header['CD2_1'] == pytest.approx(3.7219911958907E-07, abs=1e-8), f"{test_filepath} has incorrect CD2_1"
        assert hdu[hdunum].header['CD2_2'] == pytest.approx(-0.00010826103272992, abs=1e-8), f"{test_filepath} has incorrect CD2_2"
