#!/usr/bin/env bash
set -euo pipefail
test -f config/config.yaml
test -f config/samples.tsv
echo "Smoke test OK"
