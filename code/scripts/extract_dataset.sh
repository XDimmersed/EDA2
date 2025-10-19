#!/usr/bin/env bash
set -euo pipefail
repo_root="$(git rev-parse --show-toplevel)"
zip_path="${repo_root}/赛题三公开数据集.zip"
if [ ! -f "$zip_path" ]; then
  echo "[extract_dataset] 未找到压缩包：$zip_path" >&2
  exit 1
fi
if ! unzip -tqq "$zip_path" >/dev/null 2>&1; then
  cat <<MSG >&2
[extract_dataset] 压缩包无法打开，可能仍是 Git LFS 指针文件。
请先运行 'git lfs fetch' 从远端同步真实数据，再重试。
MSG
  exit 1
fi
dest_dir="$repo_root"
mkdir -p "$dest_dir"
unzip -oq "$zip_path" -d "$dest_dir"
echo "[extract_dataset] 数据集已解压到 $dest_dir"
