#!/usr/bin/env bash

# Legacy compatibility shim.
# Prefer passing --config config_local.sh or --config config_shaheen.sh.

# shellcheck source=/dev/null
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config_local.sh"
