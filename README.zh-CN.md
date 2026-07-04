# scStudio

> 🌐 **语言：** **中文** · [English README](README.md)

**在你自己的电脑上做交互式单细胞 RNA-seq 分析。** 一条命令拉起浏览器界面，上传 count 矩阵，
即可交互式完成一整套现代 scRNA-seq 流程——质控、去双细胞、归一化、降维、批次整合、聚类、
嵌入可视化、marker 检测、细胞类型注释——**全部用你自己电脑的 CPU 和内存计算，无需云服务器**。

每一步都同时照顾新手和老手：

- 💡 通俗的**「这一步是什么？」**解释卡（配一个例子）
- 🔧 **方法可选**（每步都有备选，如降维 UMAP 或 t-SNE）
- 🎚️ **可调阈值**，并给出合理默认值
- 📊 **结果总结** + **交互式预览图**（**鼠标悬停看细节**）
- 🧾 **可复现日志**，可导出成 R 脚本

> **当前状态：** 早期骨架（v0.1）。UI、模块结构、分析封装都已就位，但**尚未在装好
> Seurat/Bioconductor 的环境里端到端实跑过**——见[注意事项](#注意事项)。

---

## 为什么"localhost 优先"？

真实的 scRNA-seq 分析（Seurat/Bioconductor）需要原生计算和真实内存，纯浏览器（WASM）跑不动。
所以 scStudio 采用**本地服务 + 浏览器界面**：界面是网页，但所有计算都在**你的机器**上完成。
这和 `cellxgene launch` 是同一个模式。

---

## 三种运行方式

按你愿意安装多少东西来选。

### A. 你已经装了 R（最轻）
```r
# install.packages("remotes")
remotes::install_github("Tianqi-Ma/scStudio")
scStudio::run_app()   # 自动打开浏览器
```
需要 R 以及那些重的分析包（Seurat、Bioconductor）。下载最小，一条命令。

### B. 没有 R、零依赖 → Docker（多数用户推荐）✅
R + Seurat + Bioconductor + 应用**全部打进一个镜像**，你只需装
[Docker](https://www.docker.com/products/docker-desktop/)。
```bash
docker run --rm -p 3838:3838 -m 16g tianqima/scstudio
# 然后浏览器打开 http://localhost:3838
```
数据通过**浏览器上传**（无需挂载目录）。给 Docker 足够内存（`-m 16g`，大数据还需在
Docker Desktop 里调高内存上限）。

> 也可以自己从本仓库构建镜像：`docker build -t scstudio .`

### C. 纯小白双击即用（规划中）
用桌面安装包（Tauri/Electron/electricShine）把 R 和应用打包，用户**双击即可**——无需命令行、
无需 Docker。计划在后续版本推出。

---

## 零门槛试用（无需自己的数据）

在**导入**步骤，保持数据来源为 **Demo 数据**，点击"加载 demo 数据"即可。内置的小示例会**立即
加载（离线可用）**；你也可以选真实公共数据集（10x PBMC 1k/5k，首次使用时下载），或粘贴任意
`.rds` / `.h5` 文件的 URL 在线加载。

## 可上传的数据

- **RDS** — 保存的 `Seurat` 或 `SingleCellExperiment` 对象（旧版 Seurat 对象会自动升级）
- **10x `.h5`** — Cell Ranger 的 HDF5
- **count 表格** — CSV/TSV，基因为行、细胞为列

浏览器上传上限调得很高（默认 5 GB；用 `run_app(max_upload_mb = ...)` 修改）。

---

## 分析步骤

| # | 步骤 | 可选方法 | 你可调 |
|---|------|-----------------|-------------|
| 1 | 导入 | RDS / 10x .h5 / 表格 | 格式 |
| 2 | 质控 QC | **MAD 自适应** / 手动 | MAD 倍数或固定阈值、物种 |
| 3 | 去双细胞 | **scDblFinder** / DoubletFinder | 标记 vs 删除、分数阈值 |
| 4 | 归一化 | **LogNormalize** / SCTransform | scale factor |
| 5 | 特征选择 + PCA | HVG **vst** / mvp / dispersion | HVG 数、PC 数 |
| 6 | 批次整合 *(可选)* | **无** / Harmony / CCA / RPCA | 批次列 |
| 7 | 聚类 | **Leiden** / Louvain | 分辨率、维数 |
| 8 | 嵌入可视化 | **UMAP** / t-SNE / PaCMAP | neighbors、min_dist / perplexity |
| 9 | Marker 基因 | **wilcox** / roc / MAST | logFC、min.pct、top-N |
| 10 | 注释 | **手动** / SingleR / Azimuth | 参考集、逐簇标签 |
| 11 | 可视化 | UMAP / violin / dotplot / feature / heatmap | 基因、分组 |
| 12 | 导出 | .rds / .h5ad(尽力) / 图 / R 脚本 | — |

方法与默认值遵循当前（2023–2025）最佳实践（MAD 质控、scDblFinder、Leiden、Harmony，
pseudobulk 差异表达将在后续版本加入等）。

---

## 注意事项

- **尚未端到端验证。** 已构建并通过语法解析，但未在装好 Seurat/Bioconductor 的环境里用真实
  数据跑过。v0.1 请视为可用骨架，首次实跑时预期需要修一些细节。
- **重依赖在安装时是可选的**（放在 `Suggests`）。缺依赖时 UI 也能起；每个计算步骤会检查所需包，
  缺了会提示你安装。Docker 镜像已内置全部依赖。
- **部分步骤首次需联网**（SingleR/Azimuth 下载参考集）。
- **`.h5ad` 导出**在纯 R 下依赖 SeuratDisk（尽力而为）；主导出格式是 `.rds`。
- **内存**随细胞数增长。<10 万细胞在 16–32 GB 上较舒适；更大需更多内存。应用会警告并提供下采样。

---

## 开发

R 包结构（golem 风格）：`R/app_ui.R`、`R/app_server.R`、`R/run_app.R`，每步一个 `R/mod_*.R`，
计算封装在 `R/fct_compute.R`，通用 helper 在 `R/utils_*.R` 和 `R/fct_*.R`。开发时本地运行：
```r
pkgload::load_all(); run_app()
```

## 许可证

MIT。
