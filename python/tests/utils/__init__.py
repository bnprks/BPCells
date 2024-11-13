# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

import hashlib
import os
import shutil
import urllib.request

import pytest

# Helpers for tests that require downloading real-world data files
# Example usage:
# @conftest.slow_data_test
# def my_test(fetch_cached_file):
#   my_file = fetch_cached_file(my_url)
#   with open(my_file, "r") as f:
#      ... 

# Allow skipping if the data cache folder isn't set, but don't allow skips in a CI environment
slow_data_test = pytest.mark.skipif("BPCELLS_PYTEST_DATA_CACHE" not in os.environ and "CI" not in os.environ, reason="Data-dependent tests require BPCELLS_PYTEST_DATA_CACHE to be defined")

def fetch_cached_file(url):
    # For usage example, see: https://docs.pytest.org/en/6.2.x/fixture.html#factories-as-fixtures
    
    """Download a file from a URL and cache it for later calls.
    Cache dir is taken from os.environ["BPCELLS_PYTEST_DATA_CACHE"], and an exception is raised
    if the environment variable is not found

    Args:
        url (str): URL to fetch

    Returns:
        str: Path of downloaded file
    """
    if "BPCELLS_PYTEST_DATA_CACHE" not in os.environ:
        raise Exception("fetch_cached_file() requires environment variable BPCELLS_PYTEST_DATA_CACHE to be defined")
    
    cache_dir = os.environ["BPCELLS_PYTEST_DATA_CACHE"]
    
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)

    hash = hashlib.sha256(url.encode()).hexdigest()
    name = os.path.basename(url)
    out_path = os.path.join(cache_dir, hash + "_" + name)
    if os.path.exists(out_path):
        return out_path
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (X11; U; Linux i686) Gecko/20071127 Firefox/2.0.0.11"})
    with urllib.request.urlopen(req) as f:
        shutil.copyfileobj(f, open(out_path, "wb"))
    return out_path

