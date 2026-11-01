# 如何将项目发布到你的 GitHub

建议仓库名：`sanger-pair-assembler`。

## 方式 A：网页上传（最简单）

1. 登录 <https://github.com/new>，创建仓库 `sanger-pair-assembler`，可见性选 **Public**。
2. 不必勾选自动生成 README / .gitignore / LICENSE，因为项目中已经包含这些文件。
3. 在 Mac 上解压项目 ZIP，进入解压后的项目文件夹。
4. 在 GitHub 新仓库页面选 **uploading an existing file**，把项目根目录**内的所有文件与文件夹**上传并提交。
5. 检查 `.github/workflows/tests.yml` 是否确实上传；网页或 Finder 默认可能隐藏以 `.` 开头的目录。若无法显示隐藏目录，建议使用方式 B。
6. 在仓库页面右侧 `About` 中设置简介和 topics：`sanger-sequencing`, `16s-rrna`, `biopython`, `bioinformatics`, `consensus-sequence`, `ab1`。
7. 在 **Actions** 页面检查 CI 自动测试是否通过。第一次运行可能需要你在 GitHub 页面启用 Actions。

建议仓库简介：

> Reference-free, quality-aware forward/reverse Sanger read assembly with AB1 trimming, overlap QC, consensus FASTA and discrepancy reports.

## 方式 B：Terminal / git 上传（更稳定）

在 macOS 终端中进入项目根目录：

```bash
cd ~/Downloads/sanger-pair-assembler
git init
git add .
git commit -m "Release v0.2.0: batch assembly support"
git branch -M main
```

在 <https://github.com/new> 创建空仓库 `sanger-pair-assembler` 后，按你的 GitHub 用户名修改以下地址：

```bash
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/sanger-pair-assembler.git
git push -u origin main
```

GitHub 通常不支持通过普通登录密码执行 `git push`；可使用 GitHub CLI (`gh auth login`)、系统凭据管理器或 Personal Access Token。需要安装 GitHub CLI 时，可以使用 <https://cli.github.com/>。

## 发布 v0.2.0 Release（可选）

```bash
git tag v0.2.0
git push origin v0.2.0
```

之后在 GitHub 的 **Releases → Draft a new release** 选 `v0.2.0`，介绍功能，并附上项目 ZIP 或 wheel（如果你希望提供离线下载）。

## 发布前自查

- [ ] 仓库里没有真实 AB1、患者/个人信息、实验室内部编号和未获许可的原始序列
- [ ] `README.md`、`README.zh-CN.md`、`LICENSE`、`pyproject.toml` 均已存在
- [ ] `.github/workflows/tests.yml` 已上传，自动测试通过
- [ ] 使用模拟示例运行成功：`sanger-merge -f examples/toy_F.fasta -r examples/toy_R.fasta -o toy_merged.fasta`
- [ ] 有需要时修改许可证的署名或增加作者、机构和致谢；不要未经许可使用他人的名字或机构名

本项目以 MIT License 开源。若后续发布科研论文或作为正式研究分析工具，建议保存原始 AB1、软件版本、修剪阈值、比对参数和人工峰图复核记录。
