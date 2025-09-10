#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
覆盖分析模块：分析探针在参考基因组和FASTA文件中的覆盖情况
"""
import os
import matplotlib.pyplot as plt
import pandas as pd

# 确保中文正常显示
plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]

# 导入工具函数
from utils import read_fasta_file, read_blast_output, ensure_dir_exists

class CoverageAnalyzer:
    """探针覆盖情况分析器类"""
    
    def __init__(self, config):
        """初始化覆盖分析器"""
        self.config = config
        self.probes = {}
        self.blast_results = {}
        self.stats = {}
    
    def load_probes(self, probes_file=None):
        """加载探针序列"""
        file_path = probes_file or self.config.InputFiles.PROBES_FILE
        print(f"\n加载探针序列文件: {file_path}")
        self.probes = read_fasta_file(file_path)
        return self.probes
    
    def load_blast_results(self, blast_files=None):
        """加载BLAST输出结果"""
        files_to_load = blast_files or self.config.InputFiles.BLAST_FILES
        print(f"\n加载BLAST输出文件: {len(files_to_load)}个文件")
        
        self.blast_results = {}
        
        for file_path in files_to_load:
            print(f"  正在加载: {os.path.basename(file_path)}")
            file_results = read_blast_output(file_path)
            self.blast_results[file_path] = file_results
        
        print(f"成功加载 {len(self.blast_results)} 个BLAST文件的结果")
        return self.blast_results
    
    def analyze_coverage(self):
        """分析探针覆盖情况"""
        print(f"\n开始分析探针覆盖情况...")
        
        # 初始化统计信息
        self.stats = {
            'total_probes_in_fasta': len(self.probes),
            'total_probes_in_blast': 0,
            'probes_found_in_fasta': 0,
            'coverage_ratio': 0.0,
            'file_stats': []
        }
        
        # 收集所有blast文件中的探针ID
        all_blast_probes = set()
        
        # 分析每个blast文件
        for file_path, file_results in self.blast_results.items():
            file_name = os.path.basename(file_path)
            
            # 当前文件中的探针ID集合
            file_probes = set(file_results.keys())
            
            # 计算在FASTA文件中找到的探针数量
            found_in_fasta = file_probes.intersection(set(self.probes.keys()))
            
            # 统计信息
            file_stat = {
                'file_name': file_name,
                'total_probes': len(file_probes),
                'found_in_fasta': len(found_in_fasta),
                'not_found_in_fasta': len(file_probes) - len(found_in_fasta),
                'found_ratio': len(found_in_fasta) / len(file_probes) * 100 if file_probes else 0
            }
            
            self.stats['file_stats'].append(file_stat)
            all_blast_probes.update(file_probes)
        
        # 计算总体统计
        self.stats['total_probes_in_blast'] = len(all_blast_probes)
        self.stats['probes_found_in_fasta'] = len(all_blast_probes.intersection(set(self.probes.keys())))
        self.stats['coverage_ratio'] = self.stats['probes_found_in_fasta'] / self.stats['total_probes_in_fasta'] * 100 if self.stats['total_probes_in_fasta'] else 0
        
        # 打印分析结果
        self._print_analysis_results()
        
        return self.stats
    
    def _print_analysis_results(self):
        """打印分析结果"""
        print(f"\n===== 探针覆盖情况分析结果 =====")
        print(f"FASTA文件中的探针总数: {self.stats['total_probes_in_fasta']}")
        print(f"所有BLAST文件中的唯一探针ID总数: {self.stats['total_probes_in_blast']}")
        print(f"在FASTA文件中找到的探针总数: {self.stats['probes_found_in_fasta']}")
        print(f"在FASTA文件中未找到的探针总数: {self.stats['total_probes_in_blast'] - self.stats['probes_found_in_fasta']}")
        print(f"BLAST文件覆盖FASTA文件的比例: {self.stats['coverage_ratio']:.2f}%")
        
        print(f"\n各文件详细统计:")
        for file_stat in self.stats['file_stats']:
            print(f"  文件: {file_stat['file_name']}")
            print(f"    探针数量: {file_stat['total_probes']}")
            print(f"    在FASTA中找到: {file_stat['found_in_fasta']} ({file_stat['found_ratio']:.2f}%)")
            print(f"    在FASTA中未找到: {file_stat['not_found_in_fasta']}")
    
    def generate_report(self, output_file=None):
        """生成覆盖分析报告"""
        file_path = output_file or self.config.OutputFiles.COVERAGE_REPORT
        ensure_dir_exists(os.path.dirname(file_path))
        
        print(f"\n生成覆盖分析报告: {file_path}")
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write("===== 探针覆盖情况分析报告 =====\n\n")
            f.write(f"FASTA文件中的探针总数: {self.stats['total_probes_in_fasta']}\n")
            f.write(f"所有BLAST文件中的唯一探针ID总数: {self.stats['total_probes_in_blast']}\n")
            f.write(f"在FASTA文件中找到的探针总数: {self.stats['probes_found_in_fasta']}\n")
            f.write(f"在FASTA文件中未找到的探针总数: {self.stats['total_probes_in_blast'] - self.stats['probes_found_in_fasta']}\n")
            f.write(f"BLAST文件覆盖FASTA文件的比例: {self.stats['coverage_ratio']:.2f}%\n\n")
            
            f.write("各文件详细统计:\n")
            f.write("文件名\t探针总数\t在FASTA中找到\t在FASTA中未找到\t找到比例(%)\n")
            
            for file_stat in self.stats['file_stats']:
                f.write(
                    f"{file_stat['file_name']}\t"
                    f"{file_stat['total_probes']}\t"
                    f"{file_stat['found_in_fasta']}\t"
                    f"{file_stat['not_found_in_fasta']}\t"
                    f"{file_stat['found_ratio']:.2f}\n"
                )
        
        print(f"分析报告已保存到 {file_path}")
        return file_path
    
    def generate_visualization(self, output_file=None):
        """生成覆盖情况可视化图表"""
        file_path = output_file or self.config.OutputFiles.COVERAGE_PLOT
        ensure_dir_exists(os.path.dirname(file_path))
        
        print(f"\n生成覆盖情况可视化图表: {file_path}")
        
        # 创建两个子图
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.config.Params.FIGURE_SIZE)
        
        # 图1: 各个文件的探针分布
        file_names = [stat['file_name'] for stat in self.stats['file_stats']]
        found = [stat['found_in_fasta'] for stat in self.stats['file_stats']]
        not_found = [stat['not_found_in_fasta'] for stat in self.stats['file_stats']]
        
        ax1.bar(file_names, found, label='在FASTA中找到')
        ax1.bar(file_names, not_found, bottom=found, label='在FASTA中未找到')
        ax1.set_xlabel('BLAST输出文件')
        ax1.set_ylabel('探针数量')
        ax1.set_title('各个BLAST文件中的探针分布')
        ax1.legend()
        ax1.tick_params(axis='x', rotation=90)
        
        # 图2: 总体覆盖情况
        labels = ['FASTA文件中的探针总数', 'BLAST文件中的唯一探针总数', '在FASTA中找到的探针数']
        sizes = [
            self.stats['total_probes_in_fasta'],
            self.stats['total_probes_in_blast'],
            self.stats['probes_found_in_fasta']
        ]
        
        ax2.bar(labels, sizes, color=['blue', 'green', 'red'])
        ax2.set_ylabel('探针数量')
        ax2.set_title('探针覆盖总体情况')
        
        # 添加文本标签
        for i, v in enumerate(sizes):
            ax2.text(i, v + 100, str(v), ha='center')
        
        # 添加比例信息
        coverage_ratio = self.stats['coverage_ratio']
        ax2.text(1.5, max(sizes) * 0.8, 
                 f'BLAST覆盖FASTA的比例: {coverage_ratio:.2f}%', 
                 ha='center', 
                 bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.5))
        
        plt.tight_layout()
        plt.savefig(file_path, dpi=self.config.Params.FIGURE_DPI)
        print(f"可视化图表已保存为 {file_path}")
        
        # 关闭图表以释放内存
        plt.close(fig)
        
        return file_path
    
    def generate_score_distribution_plot(self, scores, output_file=None):
        """生成分数分布图表"""
        file_path = output_file or self.config.OutputFiles.SCORE_DISTRIBUTION_PLOT
        ensure_dir_exists(os.path.dirname(file_path))
        
        print(f"\n生成分数分布图表: {file_path}")
        
        # 从分数数据中提取分数列表
        score_values = [s['score'] for s in scores]
        
        # 创建分数分布直方图
        plt.figure(figsize=self.config.Params.FIGURE_SIZE)
        
        # 设置分数区间
        bins = [0, 20, 40, 60, 80, 100]
        plt.hist(score_values, bins=bins, edgecolor='black', alpha=0.7)
        
        # 设置图表属性
        plt.xlabel('特异性分数')
        plt.ylabel('探针数量')
        plt.title('探针特异性分数分布')
        plt.xticks(bins)
        
        # 添加网格线
        plt.grid(axis='y', alpha=0.75)
        
        # 保存图表
        plt.tight_layout()
        plt.savefig(file_path, dpi=self.config.Params.FIGURE_DPI)
        print(f"分数分布图表已保存为 {file_path}")
        
        # 关闭图表以释放内存
        plt.close()
        
        return file_path
    
    def run(self, scores=None):
        """运行完整的覆盖分析流程"""
        print("===== 探针覆盖分析流程开始 =====")
        
        # 加载数据
        self.load_probes()
        self.load_blast_results()
        
        # 分析覆盖情况
        self.analyze_coverage()
        
        # 生成报告
        self.generate_report()
        
        # 生成可视化图表
        self.generate_visualization()
        
        # 如果提供了分数数据，生成分数分布图表
        if scores:
            self.generate_score_distribution_plot(scores)
        
        print("===== 探针覆盖分析流程完成 =====")
        
        return self.stats

# 主函数（用于测试）
if __name__ == '__main__':
    # 导入配置
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from config import config, create_directories
    
    # 创建必要的目录
    create_directories()
    
    # 创建分析器实例并运行
    analyzer = CoverageAnalyzer(config)
    analyzer.run()