# Sanger Pair Assembler | Windows

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

### CMD（命令提示符）

```bat
cd /d "%USERPROFILE%\Downloads\sanger-pair-assembler"
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m sanger_assembler --input-dir "%USERPROFILE%\Desktop\SangerReads" --output-dir "%USERPROFILE%\Desktop\SangerReads\assembled_results"
```

### 只使用独立脚本

无需安装 GitHub 包，仅需安装 Biopython，然后直接运行根目录的 `sanger_merge.py`：

```powershell
.\.venv\Scripts\python.exe -m pip install biopython
.\.venv\Scripts\python.exe .\sanger_merge.py --input-dir "$HOME\Desktop\SangerReads" --output-dir "$HOME\Desktop\SangerReads\assembled_results"
```

### Conda 用户

在 Anaconda Prompt 中执行：

```bat
conda create -n sanger python=3.12 -y
conda activate sanger
cd /d "%USERPROFILE%\Downloads\sanger-pair-assembler"
python -m pip install -e .
sanger-merge --input-dir "%USERPROFILE%\Desktop\SangerReads"
```

### 常见问题

- **找不到 `py`**：改用 `python -m venv .venv`；先确认 Python 已安装。
- **PowerShell 提示禁止执行脚本**：本教程不需要激活虚拟环境，直接运行 `.\.venv\Scripts\python.exe` 即可。
- **找不到 `Bio`**：使用同一个解释器安装依赖：`.\.venv\Scripts\python.exe -m pip install biopython`。
- **路径有空格、括号或中文**：始终加双引号。不要照搬 Mac 上的 `source` 和 `~/` 命令。
- **Excel 打开 TSV 中文乱码**：在 Excel 中选择“数据 → 从文本/CSV”，按 UTF-8 编码、制表符分隔导入。
- **某些样本未成功拼接**：查看 `batch_summary.tsv`；只有成功的样本进入 `all_merged.fasta` 和 `all_merged.tsv`。

注意：自动拼接不等于人工验证。对于双向差异和插入/缺失，正式提交或发表前仍需核对原始峰图。
