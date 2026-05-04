from importlib import resources
import json
from functools import lru_cache

@lru_cache
def load_wcs_keys():
    path = resources.files("aligner").joinpath("wcs_headers.json")
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)