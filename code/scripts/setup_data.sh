#!/usr/bin/env bash
set -e

echo "Fetching LFS files (real dataset payload)..."
git lfs install
git lfs fetch --all
git lfs pull

echo "Dataset is now ready. Unpacking..."
unzip "赛题三公开数据集.zip" -d dataset
echo "Done."
