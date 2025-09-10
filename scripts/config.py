#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置文件：存储项目的路径和参数设置
"""
import os

# 项目根目录
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 数据目录
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')

# 结果目录
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')

# 报告目录
REPORTS_DIR = os.path.join(PROJECT_ROOT, 'reports')

# 输入文件路径
class InputFiles:
    # 探针序列文件
    PROBES_FILE = os.path.join(DATA_DIR, 'extracted_sequences_sliding_window.fasta')
    
    # PAV序列文件
    PAV_FILE = os.path.join(DATA_DIR, 'all_pav.fasta')
    
    # 自动检测BLAST文件和对应的基因组文件
    BLAST_FILES = []
    GENOME_FILES = []
    BLAST_GENOME_MAPPING = {}
    
    if os.path.exists(DATA_DIR):  # 确保目录存在
        # 1. 首先收集所有基因组FASTA文件
        genome_files = []
        for file_name in os.listdir(DATA_DIR):
            if file_name.endswith(('.fa', '.fasta', '.fna')) and not file_name.startswith('all_pav') and not file_name.startswith('extracted_sequences'):
                genome_path = os.path.join(DATA_DIR, file_name)
                genome_files.append((file_name, genome_path))
        
        # 2. 收集所有BLAST输出文件并建立映射
        for file_name in os.listdir(DATA_DIR):
            if file_name.startswith('extracted_blast_out_') and file_name.endswith('.txt'):
                blast_path = os.path.join(DATA_DIR, file_name)
                BLAST_FILES.append(blast_path)
                
                # 尝试从BLAST文件名提取物种标识符，例如从"extracted_blast_out_B73.txt"提取"B73"
                species_id = file_name[len('extracted_blast_out_'):-4]  # 去掉前缀和.txt后缀
                
                # 查找匹配的基因组文件
                matched_genome = None
                for genome_name, genome_path in genome_files:
                    # 如果基因组文件名包含物种标识符，则认为它们匹配
                    if species_id.lower() in genome_name.lower():
                        matched_genome = genome_path
                        if genome_path not in GENOME_FILES:
                            GENOME_FILES.append(genome_path)
                        break
                
                # 记录映射关系
                if matched_genome:
                    BLAST_GENOME_MAPPING[blast_path] = matched_genome
        
        # 3. 设置参考基因组（如果有匹配的基因组）
        if GENOME_FILES:
            REFERENCE_GENOME = GENOME_FILES[0]  # 默认使用第一个基因组作为参考
        else:
            # 如果没有找到任何基因组文件，使用默认路径
            REFERENCE_GENOME = os.path.join(DATA_DIR, 'Zm-Mo17-REFERENCE-CAU-2.0.fa')
            GENOME_FILES.append(REFERENCE_GENOME)

# 输出文件路径
class OutputFiles:
    # 评分结果文件
    SCORE_RESULTS = os.path.join(RESULTS_DIR, 'probe_scores.txt')
    
    # 高特异性探针文件
    HIGH_QUALITY_PROBES = os.path.join(RESULTS_DIR, 'high_quality_probes.fasta')
    
    # 覆盖分析结果
    COVERAGE_REPORT = os.path.join(REPORTS_DIR, 'probe_coverage_report.txt')
    
    # 探针映射结果
    MAPPING_RESULTS = os.path.join(REPORTS_DIR, 'probe_mapping_results.txt')
    MAPPING_RESULTS_CSV = os.path.join(REPORTS_DIR, 'probe_mapping_results.csv')
    MAPPING_RESULTS_HTML = os.path.join(REPORTS_DIR, 'probe_mapping_results.html')
    
    # 可视化文件
    COVERAGE_PLOT = os.path.join(REPORTS_DIR, 'probe_coverage_analysis.png')
    SCORE_DISTRIBUTION_PLOT = os.path.join(REPORTS_DIR, 'score_distribution.png')

# 分析参数
class Params:
    # 探针评分参数
    SCORE_THRESHOLD = 80  # 高特异性探针的分数阈值
    MIN_PROBE_LENGTH = 10  # 最小探针长度
    MAX_PROBE_LENGTH = 1000  # 最大探针长度
    
    # BLAST结果处理参数
    MIN_BIT_SCORE = 0  # 最小bit score，用于过滤低质量匹配
    MAX_E_VALUE = 10  # 最大e-value，用于过滤低质量匹配
    
    # 多进程参数
    MAX_WORKERS = min(16, os.cpu_count() * 2)  # 最多使用16个进程，或CPU核心数的2倍
    CHUNK_SIZE = 100  # 每批处理的探针数量
    
    # 可视化参数
    FIGURE_DPI = 300  # 图表DPI
    FIGURE_SIZE = (15, 6)  # 图表大小
    
    # 日志参数
    LOG_LEVEL = 'INFO'  # 日志级别: DEBUG, INFO, WARNING, ERROR

# 创建必要的目录
def create_directories():
    """创建项目所需的所有目录"""
    for dir_path in [DATA_DIR, RESULTS_DIR, REPORTS_DIR]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(f"创建目录: {dir_path}")

# 验证输入文件是否存在
def validate_input_files():
    """验证必要的输入文件是否存在"""
    missing_files = []
    
    # 检查必要的输入文件
    for file_path in [
        InputFiles.PROBES_FILE,
        InputFiles.PAV_FILE
    ]:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    # 检查参考基因组文件（如果指定了的话）
    if hasattr(InputFiles, 'REFERENCE_GENOME') and InputFiles.REFERENCE_GENOME and not os.path.exists(InputFiles.REFERENCE_GENOME):
        missing_files.append(InputFiles.REFERENCE_GENOME)
    
    # 检查至少有一个BLAST文件
    if not InputFiles.BLAST_FILES:
        missing_files.append("No BLAST files found in data directory")
    else:
        for blast_file in InputFiles.BLAST_FILES:
            if not os.path.exists(blast_file):
                missing_files.append(blast_file)
    
    # 检查所有指定的基因组文件是否存在
    if hasattr(InputFiles, 'GENOME_FILES'):
        for genome_file in InputFiles.GENOME_FILES:
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
    print(f"探针序列文件: {InputFiles.PROBES_FILE}")
    print(f"PAV序列文件: {InputFiles.PAV_FILE}")
    print(f"参考基因组文件: {InputFiles.REFERENCE_GENOME}")
    
    # 打印BLAST文件和基因组的映射关系
    if hasattr(InputFiles, 'BLAST_GENOME_MAPPING') and InputFiles.BLAST_GENOME_MAPPING:
        print(f"\nBLAST文件与基因组映射关系:")
        for blast_file, genome_file in InputFiles.BLAST_GENOME_MAPPING.items():
            print(f"  {os.path.basename(blast_file)} -> {os.path.basename(genome_file)}")
    
    return True

if __name__ == '__main__':
    # 如果直接运行此文件，创建目录并验证输入文件
    create_directories()
    validate_input_files()