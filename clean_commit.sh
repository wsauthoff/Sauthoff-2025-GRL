#!/bin/bash

# Install pre-commit if not already installed
pip install --quiet pre-commit

# Install the git hook (only needs to be done once per repo)
pre-commit install

# Run the hooks on all files (e.g., strip notebook outputs)
pre-commit run --all-files
