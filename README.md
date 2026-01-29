# Cytonuclear_Interaction

本仓库包含用于从 PAF 比对结果解析、分类并汇总线粒体/叶绿体插入（numt / nupt）结构证据的脚本与驱动流程。脚本以 Python 编写，配套有简单的 Bash 驱动脚本用于并行处理样本并汇总结果。

## 目录概览（与脚本相关）
- `parse_paf_and_classify.py`  
  - 主要脚本：解析 PAF 文件、基于位置映射（position mapping CSV）判定/验证候选插入区间，并为每个样本输出验证表和 QC 信息。
- `merge_sample_outputs.py`  
  - 将所有样本的 `region_validation.tsv` 合并为一个 `all_samples_region_validation.tsv`，并把各样本的 `qc.json` 合并为 `qc_all_samples.json`。
- `nupt_driver.sh` / `numt_driver.sh`  
  - 并行驱动脚本，按样本并行调用 `parse_paf_and_classify.py`，并在完成后调用 `merge_sample_outputs.py` 做最终合并。
- 其它数据目录（脚本中引用的）
  - `06_minimap2_mapping/{numt,nupt}`: 存放由 minimap2 输出的 PAF 文件，文件命名遵循 `SAMPLE_<type>_mapped.paf` 之类的模式。
  - `04_bed_by_flanking_status_v2/{numt,nupt}_position_mapping.csv`: 位置映射 CSV（由上游步骤生成），用于把 PAF 记录映射到待验证的 region metadata。
  - `07_paf_structural_validation_v1/out`: 输出目录，脚本将把 per-sample 与汇总结果写入此处。

## 功能简介
- parse_paf_and_classify.py：  
  - 对每个 PAF 文件提取 alignments（PAF record），筛选并选择最佳对齐，计算插入或重排的 “证据” 指标（覆盖率、支持的读数数、插入重叠等），为每个候选 region 产生一行 TSV 记录（`region_validation.tsv`）并写入 QC 信息 (`qc.json`)。
  - 可选生成 per-read 详情（`--write-read-details`）。
  - 支持只处理指定样本（`--samples`），支持不合并汇总（`--no-combined`）。
- merge_sample_outputs.py：  
  - 将指定类型（`numt` 或 `nupt`）下每个样本目录的 `region_validation.tsv` 合并为 `all_samples_region_validation.tsv`（保持 header 校验一致），并把所有 `qc.json` 合并为 `qc_all_samples.json`。

## 依赖
- Python 3（>=3.6）
- 使用到的 Python 库均来自标准库：`argparse`, `csv`, `glob`, `json`, `os`, `dataclasses`, `typing` 等（不需要额外 pip 包）。
- 外部：PAF 文件通常由 minimap2/其它比对工具生成；确保上游比对步骤已完成并把 PAF 文件放在指定目录。

## 使用说明

### parse_paf_and_classify.py
基本调用格式：
```
python3 parse_paf_and_classify.py \
  --type <numt|nupt> \
  --paf-dir <paf_dir> \
  --position-mapping-csv <position_mapping.csv> \
  --out-dir <out_dir> \
  [--samples "S1,S2,..."] \
  [--mapq-min N] \
  [--write-read-details] \
  [--no-combined] \
  [--include-tp P,S] \
  [--edge-window N] \
  [--min-target-cov 0.95] \
  [--min-overlap-bp N] \
  [--insert-min-overlap N] \
  [--min-support-reads N]
```

常用参数说明（默认值见括号）：
- `--type`: 必需，选择 `numt` 或 `nupt`。
- `--paf-dir`: 必需，包含 PAF 文件的目录（匹配 `*_<type>_mapped.paf`）。
- `--position-mapping-csv`: 必需，位置映射 CSV（region metadata），脚本会加载映射用于判定 region。
- `--out-dir`: 必需，输出根目录（脚本会在 `<out_dir>/<type>/<sample>/` 下写文件）。
- `--samples`: 可选，逗号分隔样本 ID（如不提供则处理 paf-dir 中所有匹配文件）。
- `--mapq-min`: 最低 MAPQ（默认 20），用于过滤低质量比对。
- `--write-read-details`: 若指定，会写出 per-read 的详细信息（会产生更多输出文件）。
- `--no-combined`: 若指定，不在脚本内写汇总文件（适用于单样本并行运行时）。
- `--include-tp`: 默认 "P"（可设置为 "P,S" 等），用于过滤 PAF 的 tp 标签。
- `--edge-window`: 默认 200（bp），边界容差值。
- `--min-target-cov`: 默认 0.95，目标序列覆盖度阈值。
- `--min-overlap-bp`: 默认 1，最小重叠碱基数。
- `--insert-min-overlap`: 默认 10000，单一 junction 情况下的插入最小重叠阈值。
- `--min-support-reads`: 默认 2，判定 region 至少需要的支持读数。

输出示例（每个 sample）：
- `<out_dir>/<type>/<sample>/region_validation.tsv` — 每行一个 region 的验证结果，包含诸如 sample, type, region_id, validated_level, extr_chrom, extr_start, extr_end, insert_length, n_full_span_reads, n_single_junction_reads, category 等字段。
- `<out_dir>/<type>/<sample>/qc.json` — 关于该样本处理的 QC / 元信息（例如 PAF 路径、统计量等）。
- 可选的 per-read 明细文件（若 `--write-read-details`）。

当在脚本内部启用合并（默认不使用 `--no-combined`），脚本会把所有收集到的 region 写入：
- `<out_dir>/<type>/all_samples_region_validation.tsv`
- `<out_dir>/<type>/qc_all_samples.json`

### merge_sample_outputs.py
基本调用格式：
```
python3 merge_sample_outputs.py --type <numt|nupt> --out-dir <out_dir> [--samples "S1,S2,..."]
```

行为：
- 在 `<out_dir>/<type>/` 下查找每个样本目录中的 `region_validation.tsv` 和 `qc.json`。
- 合并所有 `region_validation.tsv`（会校��� header 一致性）到 `<out_dir>/<type>/all_samples_region_validation.tsv`。
- 合并所有 `qc.json` 到 `<out_dir>/<type>/qc_all_samples.json`。

参数：
- `--type`: 必需，`numt` 或 `nupt`。
- `--out-dir`: 必需，输出根目录（同 parse 脚本使用）。
- `--samples`: 可选，通过逗号指定要合并的样本；不指定则合并 out-dir 下所有样本目录。

### Driver 脚本（并行处理）
仓库提供两个示例 driver：
- `nupt_driver.sh` 和 `numt_driver.sh`：示例并行处理脚本，调用 `parse_paf_and_classify.py` 对每个样本并行处理，然后调用 `merge_sample_outputs.py` 做汇总（日志会写到 `07_paf_structural_validation_v1/out/logs`）。

示例（并行度默认 4，可传递第一个参数修改）：
```
./nupt_driver.sh 4
./numt_driver.sh 4
```

脚本内部示例命令片段（来自 driver）：
```bash
python3 07_paf_structural_validation_v1/parse_paf_and_classify.py \
  --type nupt \
  --paf-dir 06_minimap2_mapping/nupt \
  --position-mapping-csv 04_bed_by_flanking_status_v2/nupt_position_mapping.csv \
  --out-dir 07_paf_structural_validation_v1/out \
  --samples "$SAMPLE" \
  --no-combined \
  > "07_paf_structural_validation_v1/out/logs/nupt_${SAMPLE}.log" 2>&1
```
并在所有样本处理完成后运行：
```bash
python3 07_paf_structural_validation_v1/merge_sample_outputs.py \
  --type nupt \
  --out-dir 07_paf_structural_validation_v1/out \
  > "$LOGDIR/nupt_MERGE.log" 2>&1
```

## 注意事项与建议
- 输入 PAF 文件命名应与脚本中 glob 模式匹配（`*_numt_mapped.paf` 或 `*_nupt_mapped.paf`），脚本会根据文件名解析样本 ID。
- position mapping CSV 的格式与字段名需满足 `load_position_mapping` 的解析逻辑（请参考 `parse_paf_and_classify.py` 中的 `load_position_mapping` 函数以确保列和 region_id 的格式一致）。
- 并行运行时，建议在 driver 脚本中使用 `--no-combined`，待所有并行任务完成后再调用 `merge_sample_outputs.py` 做合并（驱动脚本示例已体现此做法）。
- 如果希望调试单个样本，可直接用 `parse_paf_and_classify.py` 的 `--samples` 指定单个样本并去掉 `--no-combined` 以生成合并文件（适用于仅有少量样本的场景）。

## 输出字段说明（简要）
`region_validation.tsv` 包含（按脚本中定义的顺序，示例）：
- sample, type, region_id, validated_level, extr_chrom, extr_start, extr_end, extr_length, left_flank, right_flank, insert_length, tlen, n_full_span_reads, n_single_junction_reads, max_target_cov, max_insert_overlap, n_read_region_pairs, orig_chrom, orig_start, orig_end, orig_length, category
（具体字段可在 `parse_paf_and_classify.py` 的 `region_header` 中查看并保持一致）

`qc.json` 含有 per-sample 的 QC 字段（如 paf_path、统计量等），`qc_all_samples.json` 是多个 `qc.json` 的列表。

