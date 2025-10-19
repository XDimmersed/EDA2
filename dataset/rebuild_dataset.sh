#!/bin/bash
# 合并所有分卷为完整 ZIP
cat dataset_split.z* dataset_split.zip > full_dataset.zip
# 解压完整数据
unzip full_dataset.zip
