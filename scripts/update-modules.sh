#!/bin/bash
# Helper script to update the module symlink after fmp builds
# This ensures VS Code language server can always find the latest modules

cd "$(dirname "$0")/.."

# Find the current gfortran module directory
GFORTRAN_DIR=$(find build -maxdepth 1 -name "gfortran_*" -type d | head -1)

if [ -z "$GFORTRAN_DIR" ]; then
    echo "No gfortran module directory found in build/"
    exit 1
fi

GFORTRAN_BASENAME=$(basename "$GFORTRAN_DIR")

# Update the symlink
echo "Updating build/modules -> $GFORTRAN_BASENAME"
ln -sf "$GFORTRAN_BASENAME" build/modules

echo "Module symlink updated successfully!"
echo "Modules available at: build/modules/"
ls -la build/modules/*.mod 2>/dev/null | head -5