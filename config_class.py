#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置类定义：定义项目中使用的配置对象
"""
import os

# 项目根目录（根据config.py的位置确定）
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# 输入文件路径类 - 作为顶级类，而不是嵌套类
class InputFiles:
    def __init__(self, data_dir):
        # 探针序列文件
        self.PROBES_FILE = os.path.join(data_dir, 'extracted_sequences_sliding_window_clean.fasta')
        
        # PAV序列文件
        self.PAV_FILE = os.path.join(data_dir, 'all_pav.fasta')
        
        # 自动检测BLAST文件和对应的基因组文件
        self.BLAST_FILES = []
        self.GENOME_FILES = []
        self.BLAST_GENOME_MAPPING = {}
        
        if os.path.exists(data_dir):  # 确保目录存在
            # 1. 首先收集所有基因组FASTA文件
            genome_files = []
            for file_name in os.listdir(data_dir):
                if file_name.endswith(('.fa', '.fasta', '.fna')) and not file_name.startswith('all_pav') and not file_name.startswith('extracted_sequences'):
                    genome_path = os.path.join(data_dir, file_name)
                    genome_files.append((file_name, genome_path))
            
            # 2. 收集所有BLAST输出文件并建立映射
            for file_name in os.listdir(data_dir):
                if file_name.startswith('extracted_blast_out_') and file_name.endswith('.txt'):
                    blast_path = os.path.join(data_dir, file_name)
                    self.BLAST_FILES.append(blast_path)
                    
                    # 尝试从BLAST文件名提取物种标识符，例如从"extracted_blast_out_B73.txt"提取"B73"
                    species_id = file_name[len('extracted_blast_out_'):-4]  # 去掉前缀和.txt后缀
                    
                    # 查找匹配的基因组文件
                    matched_genome = None
                    for genome_name, genome_path in genome_files:
                        # 如果基因组文件名包含物种标识符，则认为它们匹配
                        if species_id.lower() in genome_name.lower():
                            matched_genome = genome_path
                            if genome_path not in self.GENOME_FILES:
                                self.GENOME_FILES.append(genome_path)
                            break
                    
                    # 记录映射关系
                    if matched_genome:
                        self.BLAST_GENOME_MAPPING[blast_path] = matched_genome
                    else:
                        # 如果没有找到匹配的基因组，使用默认基因组（如果有）
                        print(f"警告: 未找到与BLAST文件 {file_name} 匹配的基因组文件")
            
            # 3. 设置参考基因组（如果有匹配的基因组）
            if self.GENOME_FILES:
                self.REFERENCE_GENOME = self.GENOME_FILES[0]  # 默认使用第一个基因组作为参考
            else:
                # 如果没有找到任何基因组文件，使用默认路径
                self.REFERENCE_GENOME = os.path.join(data_dir, 'Zm-Mo17-REFERENCE-CAU-2.0.fa')
                self.GENOME_FILES.append(self.REFERENCE_GENOME)

# 输出文件路径类 - 作为顶级类，而不是嵌套类
class OutputFiles:
    def __init__(self, results_dir, reports_dir):
        # 评分结果文件
        self.SCORE_RESULTS = os.path.join(results_dir, 'probe_scores.txt')
        
        # 高特异性探针文件
        self.HIGH_QUALITY_PROBES = os.path.join(results_dir, 'high_quality_probes.fasta')
        
        # 覆盖分析结果
        self.COVERAGE_REPORT = os.path.join(reports_dir, 'probe_coverage_report.txt')
        
        # 探针映射结果
        self.MAPPING_RESULTS = os.path.join(reports_dir, 'probe_mapping_results.txt')
        self.MAPPING_RESULTS_CSV = os.path.join(reports_dir, 'probe_mapping_results.csv')
        self.MAPPING_RESULTS_HTML = os.path.join(reports_dir, 'probe_mapping_results.html')
        
        # 可视化文件
        self.COVERAGE_PLOT = os.path.join(reports_dir, 'probe_coverage_analysis.png')
        self.SCORE_DISTRIBUTION_PLOT = os.path.join(reports_dir, 'score_distribution.png')

# 分析参数类 - 作为顶级类，而不是嵌套类
class Params:
    def __init__(self):
        # 探针评分参数
        self.SCORE_THRESHOLD = 80  # 高特异性探针的分数阈值
        self.MIN_PROBE_LENGTH = 10  # 最小探针长度
        self.MAX_PROBE_LENGTH = 1000  # 最大探针长度
        
        # BLAST结果处理参数
        self.MIN_BIT_SCORE = 0  # 最小bit score，用于过滤低质量匹配
        self.MAX_E_VALUE = 10  # 最大e-value，用于过滤低质量匹配
        
        # 多进程参数
        self.MAX_WORKERS = min(16, os.cpu_count() * 2)  # 最多使用16个进程，或CPU核心数的2倍
        self.CHUNK_SIZE = 100  # 每批处理的探针数量
        
        # 可视化参数
        self.FIGURE_DPI = 300  # 图表DPI
        self.FIGURE_SIZE = (15, 6)  # 图表大小
        
        # 日志参数
        self.LOG_LEVEL = 'INFO'  # 日志级别: DEBUG, INFO, WARNING, ERROR

# 定义配置类
class Config:
    """项目配置类"""
    
    def __init__(self):
        """初始化配置"""
        # 数据目录
        self.DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
        
        # 结果目录
        self.RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
        
        # 报告目录
        self.REPORTS_DIR = os.path.join(PROJECT_ROOT, 'reports')
        
        # 初始化外部定义的类实例
        self.InputFiles = InputFiles(self.DATA_DIR)
        self.OutputFiles = OutputFiles(self.RESULTS_DIR, self.REPORTS_DIR)
        self.Params = Params()

# 创建全局配置对象
config = Config()

# 创建必要的目录
def create_directories():
    """创建项目所需的所有目录"""
    for dir_path in [config.DATA_DIR, config.RESULTS_DIR, config.REPORTS_DIR]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(f"创建目录: {dir_path}")

# 验证输入文件是否存在
def validate_input_files():
    """验证必要的输入文件是否存在"""
    missing_files = []
    
    # 检查必要的输入文件
    for file_path in [
        config.InputFiles.PROBES_FILE,
        config.InputFiles.PAV_FILE
    ]:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    # 检查参考基因组文件（如果指定了的话）
    if hasattr(config.InputFiles, 'REFERENCE_GENOME') and config.InputFiles.REFERENCE_GENOME and not os.path.exists(config.InputFiles.REFERENCE_GENOME):
        missing_files.append(config.InputFiles.REFERENCE_GENOME)
    
    # 检查至少有一个BLAST文件
    if not config.InputFiles.BLAST_FILES:
        missing_files.append("No BLAST files found in data directory")
    
    # 检查所有指定的基因组文件是否存在
    if hasattr(config.InputFiles, 'GENOME_FILES'):
        for genome_file in config.InputFiles.GENOME_FILES:
            if not os.path.exists(genome_file):
                missing_files.append(genome_file)
    
    # 如果有缺失的文件，打印警告
    if missing_files:
        print("警告: 以下必要的输入文件缺失:")
        for file_path in missing_files:
            print(f"  {file_path}")
        print("请确保在运行分析前准备好这些文件")
        return False
    
    # 打印检测到的BLAST和基因组文件信息
    print(f"\n检测到的文件信息:")
    print(f"探针序列文件: {config.InputFiles.PROBES_FILE}")
    print(f"PAV序列文件: {config.InputFiles.PAV_FILE}")
    print(f"参考基因组文件: {config.InputFiles.REFERENCE_GENOME}")
    
    # 打印BLAST文件和基因组的映射关系
    if hasattr(config.InputFiles, 'BLAST_GENOME_MAPPING') and config.InputFiles.BLAST_GENOME_MAPPING:
        print(f"\nBLAST文件与基因组映射关系:")
        for blast_file, genome_file in config.InputFiles.BLAST_GENOME_MAPPING.items():
            print(f"  {os.path.basename(blast_file)} -> {os.path.basename(genome_file)}")
    
    return True

# 主函数（用于测试）
if __name__ == '__main__':
    # 打印配置信息
    print("项目配置信息:")
    print(f"数据目录: {config.DATA_DIR}")
    print(f"结果目录: {config.RESULTS_DIR}")
    print(f"报告目录: {config.REPORTS_DIR}")
    print(f"\n输入文件:")
    print(f"探针序列文件: {config.InputFiles.PROBES_FILE}")
    print(f"参考基因组文件: {config.InputFiles.REFERENCE_GENOME}")
    print(f"PAV序列文件: {config.InputFiles.PAV_FILE}")
    print(f"BLAST文件数量: {len(config.InputFiles.BLAST_FILES)}")
    print(f"\n参数设置:")
    print(f"高特异性阈值: {config.Params.SCORE_THRESHOLD}")
    print(f"最大进程数: {config.Params.MAX_WORKERS}")
    
    # 创建目录并验证文件
    create_directories()
    validate_input_files()