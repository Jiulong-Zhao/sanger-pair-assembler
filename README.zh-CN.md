# Sanger Pair Assembler（Sanger 一代测序双向拼接工具）

中文 | [English](README.md)

**本项目使用 Python/Biopython，将 F 端和 R 端的一代测序（Sanger）结果进行质量修剪、反向互补和重叠拼接，生成共识序列 FASTA。** 无需参考基因组，适合 16S rRNA、ITS、PCR 扩增子等有重叠区域的双向测序。

> 项目目前为 v0.2.1：自动共识并不等于已经人工验证。特别是测序峰图中的插入/缺失、混合峰、重复序列和引物序列，仍需核对。**不会自动去除 PCR 引物，也不适合无重叠读数的强行连接。**

## 功能

- 输入：`.ab1`、`.abi`、`.fasta`、`.fa`、`.fastq`、`.fq`、`.seq`、`.txt` 等。
- 自动对**原始方向的 R 端**反向互补，再与 F 端进行半全局重叠比对。
- 对 AB1 采用与 Biopython `abi-trim` 默认行为一致的 Modified Mott 末端修剪算法，并通过 `--trim-cutoff` 自定义阈值。
- 报告 F/R 的原始长度、质量修剪后的长度，以及原始读数 **5′ / 3′ 分别截去多少 bp**。
- 允许少量碱基冲突及 indel；有 Phred 质量分数时用于冲突判定，否则保留 IUPAC 模糊碱基；所有冲突均输出 TSV 供峰图复核。
- 无可信重叠时停止拼接，不虚构中间序列。
- 单样本输出 FASTA、F 和反向互补 R、比对文本、差异 TSV、JSON 质控报告。
- **批量模式**：自动识别文件夹内同一样本的 F/R，保留每个样本的完整输出，并额外输出汇总 FASTA、两列 TSV 与批次质控报告。

## 一、macOS 安装（含 Apple Silicon）

本项目支持 Python 3.8+，并依赖 Biopython 1.83 或以上；**Python 3.8 已停止官方安全维护**，新安装推荐用更新的 Python。若你已有名为 `python3` 的 Conda 环境（即终端前缀 `(python3)`），可以继续使用：

```bash
conda activate python3
cd ~/Downloads/sanger-pair-assembler
python -m pip install -e .
python -c "from Bio import SeqIO; print('Biopython OK')"
```

也可新建虚拟环境：

```bash
cd ~/Downloads/sanger-pair-assembler
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

**重要：** macOS 上 `python` 和 `python3` 可能不属于同一个解释器。始终使用 `python -m pip` 安装，并使用同一个 `python` 运行。

## Windows 10/11 安装和使用

Windows 与 macOS 使用同一套 Python 拼接算法，支持 AB1、FASTA、单样本和批量模式。建议使用 Python 3.12 或 3.13。下载 [Python](https://www.python.org/downloads/windows/) 并安装后，打开 PowerShell，先检查 `py --version`；如果没有 `py` 命令，可以试试 `python --version`。

### PowerShell（推荐）

下载 GitHub 项目 ZIP 并解压到“下载”文件夹，执行：

```powershell
cd "$HOME\Downloads\sanger-pair-assembler"
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m sanger_assembler --help
```

**批量拼接**：假设 AB1 都在桌面的 `SangerReads` 文件夹中：

```powershell
.\.venv\Scripts\python.exe -m sanger_assembler --input-dir "$HOME\Desktop\SangerReads" --output-dir "$HOME\Desktop\SangerReads\assembled_results" --trim-cutoff 0.05
```

**单个样本**：

```powershell
.\.venv\Scripts\python.exe -m sanger_assembler -f "$HOME\Desktop\SangerReads\S1_F.ab1" -r "$HOME\Desktop\SangerReads\S1_R.ab1" -o "$HOME\Desktop\SangerReads\S1_merged.fasta"
```

[Windows 完整使用指南](WINDOWS_GUIDE.zh-CN.md)

---

## 二、如何拼接

假设 F、R 文件在当前工作目录：

```bash
sanger-merge \
  -f 'sample_F_27F.ab1' \
  -r 'sample_R_1492R.ab1' \
  -o S1_merged.fasta
```

也可以这样运行：

```bash
python -m sanger_assembler -f 'F.ab1' -r 'R.ab1' -o result.fasta
```

**R 端要提供测序公司原始交付的文件，不要事先反向互补。** 你的文件名中可能有括号、方括号等符号，建议用单引号引住文件名。

仓库内置不含任何真实样本信息的**人工模拟示例**：

```bash
sanger-merge -f examples/toy_F.fasta -r examples/toy_R.fasta -o toy_merged.fasta
```

## 三、批量拼接整个文件夹

最适合你目前的 Mac 使用习惯：**只下载仓库根目录的 `sanger_merge.py` 即可**，用下面两行替换 `~/software` 里的旧脚本，然后对测序文件夹运行：

```bash
cp ~/Downloads/sanger_merge.py ~/software/sanger_merge.py
cd '/你的/40#测序结果'
python ~/software/sanger_merge.py --input-dir . --output-dir ./assembled_results
```

如果已通过 `python -m pip install -e .` 安装 GitHub 项目，也可以直接执行：

```bash
sanger-merge --input-dir '/path/to/ab1_folder' --output-dir '/path/to/results'
```

质量修剪仍可调整，例如：

```bash
python ~/software/sanger_merge.py -i . --output-dir ./q20_results --trim-cutoff 0.01
```

默认**仅扫描指定文件夹当前层**；若每个样本在子文件夹中，增加 `--recursive`。如果重复运行同一个输出目录，必须使用 `--force`，以免覆盖上一次的结果。

### 文件自动配对规则

支持以下文件命名方式，文件名后缀可为 `.ab1`、`.abi`、`.fasta`、`.fastq`、`.seq`、`.txt` 等：

| F 文件 | R 文件 | 自动识别样本 |
|---|---|---|
| `S1_F.ab1` | `S1_R.ab1` | `S1` |
| `S2_27F.ab1` | `S2_1492R.ab1` | `S2` |
| `sample_F_27F.ab1` | `sample_R_1492R.ab1` | `sample` |
| `S3_forward.fasta` | `S3_reverse.fasta` | `S3` |
| `0001_316_(S1)_[27F](1).ab1` | `0002_316_(S1)_[1492R](1).ab1` | `S1` |

区分大小写不影响 F/R 后缀识别。原则是**后缀必须明确标记 F/R**：`_F`、`_R`、`_27F`、`_1492R`、`_forward`、`_reverse`，或方括号、括号包裹这些标记。若文件包含 `(S1)` 形式的样本标记，优先提取 `S1`，因此能忽略测序公司不同的 `0001_` / `0002_` 文件序号。

- 缺少 F 或 R：标记为 `unpaired`，不拼接。
- 有重复的 F 或 R 文件：标记为 `ambiguous`，不随意挑选。
- 文件名无法识别：记录到 `unrecognized_files.tsv`，不会把不相关文件误认为样本。
- 即使 F/R 配对成功，重叠区 QC 不合格仍不写入总 FASTA/TSV。

### 批量输出结构

```text
assembled_results/
├── all_merged.fasta               # 所有成功拼接的序列，按样本命名
├── all_merged.tsv                 # 恰好两列：样本名称、拼接后序列
├── batch_summary.tsv              # 每个样本的成功/失败原因、修剪长度和重叠信息
├── unrecognized_files.tsv         # 未能识别 F/R 后缀的文件
└── samples/
    ├── S1/
    │   ├── S1_merged.fasta
    │   ├── S1_merged.F_trimmed.fasta
    │   ├── S1_merged.R_revcomp.fasta
    │   ├── S1_merged.alignment.txt
    │   ├── S1_merged.disagreements.tsv
    │   └── S1_merged.summary.json
    └── S2/
        └── ...
```

`all_merged.tsv` **严格只有两列**，编码为 UTF-8，每个样本一行，第二列为不换行的完整序列：

```tsv
样本名称	拼接后序列
S1	ACTG...
S2	GGCA...
```

`all_merged.fasta` 的标题为安全化的样本名称（例如 `>S1`）；TSV 则保留原本识别出的样本名称。汇总文件只包含成功拼接的样本。若有失败或缺失配对，仍生成汇总文件，但程序退出码为 1；一个成功样本都没有则退出码为 2。

**注意：** 不同目录中相同的样本名会被识别为重复文件并拒绝自动配对。确实需要拼接这类特殊命名时，建议先重命名为明确且唯一的 `样本_F`、`样本_R`。

---

## 四、质量修剪参数怎么选？

`--trim-cutoff` 表示 Modified Mott 累积质量评分中的**错误概率参数**，而非「逐个删除低于对应 Q 值的碱基」。默认 `0.05` 与 Biopython `abi-trim` 的修剪切片逻辑一致。

| 参数 | 对应单碱基 Q 值分界 | 说明 |
|---|---|---|
| `--trim-cutoff 0.05` | 约 Q13 | 默认，更容易保留双向重叠信息 |
| `--trim-cutoff 0.01` | Q20 | 更严格，建议用于敏感性比较 |
| `--trim-cutoff 0.001` | Q30 | 很严格，可能剪掉全部重叠区 |

**不建议盲目统一使用 Q30。** 对 16S 双向数据，建议先用默认 `0.05` 获得初步共识，再试 `0.01` 查看结果是否稳定，并结合原始峰图复核差异位点。Q20 通常是碱基质量评价中的常用参考，但并没有适用于所有 Sanger 项目的强制统一末端剪切阈值。

```bash
# 默认（Q13 左右的 Mott 得分阈值）
sanger-merge -f F.ab1 -r R.ab1 -o result_q13.fasta --trim-cutoff 0.05

# 更严格（Q20 得分阈值）
sanger-merge -f F.ab1 -r R.ab1 -o result_q20.fasta --trim-cutoff 0.01

# 不执行 AB1 自动质量修剪
sanger-merge -f F.ab1 -r R.ab1 -o raw_merge.fasta --no-trim
```

**说明：** 默认只对 AB1/ABI 做末端质量修剪；FASTA 本身不含质量值。FASTQ 会读取 Phred 质量值并用于重叠冲突判断，但暂不执行自动末端质量修剪。`--min-q 20` 属于冲突判定阈值，不是末端修剪阈值。

## 五、单样本结果文件

以 `-o S1_merged.fasta` 为例：

| 输出 | 内容 |
|---|---|
| `S1_merged.fasta` | 最终初步共识序列 |
| `S1_merged.F_trimmed.fasta` | 剪切后的 F 序列 |
| `S1_merged.R_revcomp.fasta` | 剪切并反向互补后的 R 序列 |
| `S1_merged.alignment.txt` | 重叠区域序列比对 |
| `S1_merged.disagreements.tsv` | 碱基冲突及插入/缺失清单、原始测序位置、Phred 质量值 |
| `S1_merged.summary.json` | 各端修剪长度、阈值、重叠区、一致性、拼接状态 |

报告中的位置从 **1** 开始计数。`R_raw_pos_1based` 是测序公司交付的原始 R 读数坐标；`R_revcomp_pos_1based` 是反向互补后的位置。单侧的插入/缺失不能仅凭另一端不存在就认定真假，应检查峰图。

为了避免覆盖分析结果，重复运行同名输出时默认报错。要明确覆盖，请加 `--force`。如果重叠区不满足阈值，则保存独立读数和失败原因，不生成虚假的拼接 FASTA。

## 六、更多命令行选项

```text
--min-overlap 40         重叠区至少 40 个比对列
--min-identity 0.90      重叠区精确匹配比例至少 90%
--min-q 20              冲突碱基判定的最低质量分数
--quality-delta 5       两方向质量差至少 5 才优先采纳高质量碱基
--input-dir DIR         从文件夹批量读取 F/R 测序文件
--output-dir DIR        批量输出目录（默认 input_dir/assembled_results）
--recursive             递归读取子文件夹
--no-trim               不修剪 AB1 末端
--force                 覆盖已有输出
```

完整帮助：`sanger-merge --help`。

## 七、测试、引用和许可

运行测试：

```bash
python -m pip install -e . pytest
python -m pytest -q
```

公开仓库仅包含**合成示例数据**，不包含真实的 S1 AB1、样本序列或个人路径。请勿在公开 GitHub Issue 上传未经授权的临床、患者、未发表或其它敏感测序数据。

本项目基于 MIT License 开源，可自由使用、修改与二次开发，但使用者需对最终生物学结论进行人工核验。

**参考资料**：

- Crossley BM et al. (2020). [Guidelines for Sanger sequencing and molecular assay monitoring](https://doi.org/10.1177/1040638720905833).
- Ewing B, Green P (1998). [Base-calling of automated sequencer traces using phred. II. Error probabilities](https://doi.org/10.1101/gr.8.3.186).
- [Biopython 官方文档](https://biopython.org/docs/latest/).
