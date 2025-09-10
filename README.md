# 玉米探针分析项目

## 项目概述
这是一个重新设计的玉米探针分析项目，专注于探针特异性评估、覆盖分析和探针映射功能。该项目提供了一套完整的工具集，用于处理、分析和可视化玉米基因芯片探针数据。

## 项目结构

```
玉米_probe_analysis/
├── scripts/           # 主脚本和功能模块
├── data/              # 输入数据目录
├── results/           # 结果文件目录
├── reports/           # 报告和可视化结果目录
├── README.md          # 项目说明文档
├── config_class.py    # 配置类定义
├── config.py          # 配置文件
├── start_analysis.py  # Python启动脚本
└── 故障排除指南.md    # 问题解决指南
```

## 目录说明

- **scripts/**: 包含所有Python脚本和功能模块
  - `main.py`: 主入口脚本，协调整个分析流程
  - `probe_scoring.py`: 探针评分模块（使用改进算法）
  - `coverage_analysis.py`: 覆盖分析模块
  - `probe_mapping.py`: 探针映射模块
  - `utils.py`: 通用工具函数

- **data/**: 存放输入数据文件
  - 探针序列FASTA文件
  - BLAST输出文件
  - 参考基因组FASTA文件
  - PAV序列FASTA文件

- **results/**: 存储分析结果文件
  - 探针评分结果
  - 高特异性探针列表

- **reports/**: 存储生成的报告和可视化结果
  - 文本报告
  - 图表文件
  - HTML报告

## 安装和环境配置

### 必要依赖
- Python 3.6+ 
- 主要Python库：
  - Biopython
  - matplotlib
  - pandas
  - numpy
  - chardet

### 安装方法
```bash
pip install biopython matplotlib pandas numpy chardet
```

## 使用方法

### 1. 准备数据
将所有需要的输入文件放入`data/`目录中：
- 探针序列文件：`extracted_sequences_sliding_window.fasta`
- 参考基因组文件：`Zm-Mo17-REFERENCE-CAU-2.0.fa`
- PAV序列文件：`all_pav.fasta`
- BLAST输出文件：`extracted_blast_out_*.txt`

### 2. 运行程序（推荐方法）
直接双击运行项目根目录下的`start_analysis.py`文件。

### 3. 命令行运行方式
```bash
cd e:\laso\玉米_probe_analysis
python start_analysis.py
```

## 分析流程

程序会自动执行以下完整的分析流程：

1. **探针评分**：使用基于非特异性指数的改进算法计算每个探针的特异性分数
2. **覆盖分析**：分析探针在参考基因组中的覆盖情况，生成统计报告和可视化图表
3. **探针映射**：将高特异性探针映射回PAV区域，生成多种格式的报告

## 输出结果

程序执行完成后，会在以下目录生成结果文件：

- **results/** 目录：
  - `probe_scores.txt`: 探针评分结果
  - `high_quality_probes.fasta`: 高特异性探针序列

- **reports/** 目录：
  - `probe_coverage_report.txt`: 探针覆盖情况报告
  - `probe_mapping_results.txt`: 探针映射文本报告
  - `probe_mapping_results.csv`: 探针映射CSV报告（便于Excel处理）
  - `probe_mapping_results.html`: 探针映射HTML报告（可视化界面）
  - `probe_coverage_analysis.png`: 探针覆盖情况图表
  - `score_distribution.png`: 探针评分分布图表

## 评分算法说明

本项目使用了改进的探针特异性评分算法，基于以下步骤：

1. **计算原始非特异性指数**：Raw_Non_Specificity = Σ(nident_i / L)
2. **计算原始特异性分数**：Raw_Specificity_Score = 1 / (1 + Raw_Non_Specificity)
3. **标准化与最终分数**：Final_Specificity_Score = 100 * Raw_Specificity_Score

这种算法能够更准确地评估探针的特异性，避免因只关注完全匹配而忽略其他高相似度匹配的问题。

## 注意事项

- 确保所有输入文件的格式正确
- 大型数据集可能需要较长的处理时间
- 可以根据需要调整`config_class.py`中的参数来优化分析结果
- 如遇到运行问题，请参考`故障排除指南.md`文件

## 许可证

本项目仅供学术研究使用。