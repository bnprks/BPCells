## BPCells Development guide

BPCells has switched to a multi-language repository setup, meaning that R and Python packages
are saved in the same repository. For the most part, single-language changes just go in their
respective repositories, though adding C++ files requires a few updates (described below)

Helper scripts included for building docs and running tests. Within `r`, `python`, or `cpp`
run `bash scripts/run_tests.sh` or `bash scripts/build_docs.sh`. Run from the root directory
to run tests/docs for all languages in one command.

## Building docs
Warning: Docs setup may change as python and C++ docs get created

1. Docs site is built into the root folder of the docs-html branch
2. Run `git worktree add r/docs docs-html` to put a copy of the docs-html branch in the r/docs folder
3. In `r` folder, run `pkgdown::build_site()`
4. In `python` folder, run `sphinx-build -M html docs/source docs/build`
    - Optionally run `rm -r docs/build/html` to get a clean re-generation
5. Run `cp -r docs/build/html/* ../r/docs/python`
6. Pull up website in local browser to check it's okay
7. In `r/docs` folder, run git commit

## R


**Dev requirements**: `install.packages(c("pkgdown", "devtools", "testthat"))`
 - Documentation uses `pkgdown`
 - Tests use `testthat`

**Installation**: `remotes::install_github("bnprks/BPCells/r")`

## Python

**Dev requirements**: [tox](https://tox.wiki/en/stable/installation.html)
 - Documentation uses `Sphinx` with [PyData theme](https://pydata-sphinx-theme.readthedocs.io/en/latest/)
 - Tests use `pytest`
 - Includes a helper file with the overal project home page, linking to C++, R, and Python docs.
 - `tox` handles virtualenv setup for tests + documentation

**Installation**: `pip install "git+https://github.com/bnprks/BPCells#subdirectory=python"`

## C++

**Dev requirements**: cmake, doxygen
 - Documentation uses Doxygen with the [Doxygen Awesome theme](https://jothepro.github.io/doxygen-awesome-css/index.html)
 - Tests use googletest + cmake + ctest

**Adding files**:
 - If C++ tests are added: Update `cpp/CMakeLists.txt`
 - If a new `.cpp` file is added: Update `r/src/Makevars.in`

### VS Code editing extensions

- Python: autoDocstring
- C++: Doxygen Documentation Generator

To avoid multiple copies of results for C++ files and searches, add to `.vscode/settings.json`:
```json
"search.exclude": {
    "r/src/bpcells-cpp/**": true,
    "python/src/bpcells-cpp/**": true,
    "r/src/vendor/**": true,
    "python/src/vendor/**": true
}
```
