# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

import pytest

import utils

@pytest.fixture
def fetch_cached_file():
    # For usage example, see: https://docs.pytest.org/en/6.2.x/fixture.html#factories-as-fixtures
    return utils.fetch_cached_file
