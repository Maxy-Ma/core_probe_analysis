#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
玉米探针分析项目主脚本

该脚本协调整个探针分析流程，包括：
1. 探针评分：使用改进算法计算探针特异性分数
2. 覆盖分析：分析探针在参考基因组中的覆盖情况
3. 探针映射：将高特异性探针映射回PAV区域

用法：
    python main.py [--steps STEPS] [--config CONFIG]

参数：
    --steps: 指定要运行的步骤，可选值为'score', 'coverage', 'mapping'，多个步骤用逗号分隔，默认为运行所有步骤
    --config: 自定义配置文件路径
"""
import os
import sys
import argparse
import time
from datetime import datetime

# 添加当前目录到Python路径
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# 导入配置和工具函数
from config import config, create_directories, validate_input_files
from utils import setup_logging, format_run_time

# 导入各个模块
from probe_scoring import ProbeScorer
from coverage_analysis import CoverageAnalyzer
from probe_mapping import ProbeMapper

class ProbeAnalysisPipeline:
    """探针分析流程管理类"""
    
    def __init__(self, steps=None, custom_config=None):
        """初始化分析流程"""
        # 使用自定义配置或默认配置
        self.config = custom_config or config
        
        # 设置步骤
        default_steps = ['score', 'coverage', 'mapping']
        self.steps = steps or default_steps
        
        # 验证步骤
        valid_steps = {'score', 'coverage', 'mapping'}
        for step in self.steps:
            if step not in valid_steps:
                raise ValueError(f"无效的步骤: {step}，有效步骤为: {', '.join(valid_steps)}")
        
        # 设置日志
        self.logger = setup_logging(self.config.Params.LOG_LEVEL)
        
        # 初始化各个模块
        self.scorer = None
        self.analyzer = None
        self.mapper = None
        
        # 存储结果
        self.scores = None
        
    def initialize_pipeline(self):
        """初始化分析流程"""
        print("===== 玉米探针分析流程初始化 ======")
        
        # 创建必要的目录
        create_directories()
        
        # 验证输入文件
        print("\n验证输入文件...")
        if not validate_input_files():
            print("警告: 部分输入文件缺失，分析可能无法正常进行")
        
        # 初始化各个模块
        if 'score' in self.steps or 'coverage' in self.steps or 'mapping' in self.steps:
            self.scorer = ProbeScorer(self.config)
        
        if 'coverage' in self.steps:
            self.analyzer = CoverageAnalyzer(self.config)
        
        if 'mapping' in self.steps:
            self.mapper = ProbeMapper(self.config)
        
        print("\n分析流程初始化完成！")
        print(f"将执行的步骤: {', '.join(self.steps)}")
    
    def run_probe_scoring(self):
        """运行探针评分步骤"""
        if 'score' not in self.steps:
            return None
        
        print("\n" + "="*50)
        print("开始探针评分步骤")
        print("="*50)
        
        start_time = time.time()
        
        # 运行评分流程
        self.scores = self.scorer.run()
        
        elapsed_time = time.time() - start_time
        print(f"\n探针评分步骤完成！")
        print(f"总耗时: {format_run_time(elapsed_time)}")
        
        return self.scores
    
    def run_coverage_analysis(self):
        """运行覆盖分析步骤"""
        if 'coverage' not in self.steps:
            return None
        
        print("\n" + "="*50)
        print("开始覆盖分析步骤")
        print("="*50)
        
        start_time = time.time()
        
        # 运行覆盖分析流程
        stats = self.analyzer.run(scores=self.scores)
        
        elapsed_time = time.time() - start_time
        print(f"\n覆盖分析步骤完成！")
        print(f"总耗时: {format_run_time(elapsed_time)}")
        
        return stats
    
    def run_probe_mapping(self):
        """运行探针映射步骤"""
        if 'mapping' not in self.steps:
            return None
        
        print("\n" + "="*50)
        print("开始探针映射步骤")
        print("="*50)
        
        start_time = time.time()
        
        # 运行映射流程
        mapping_results = self.mapper.run()
        
        elapsed_time = time.time() - start_time
        print(f"\n探针映射步骤完成！")
        print(f"总耗时: {format_run_time(elapsed_time)}")
        
        return mapping_results
    
    def run(self):
        """运行完整的分析流程"""
        print(f"\n===== 玉米探针分析流程开始 (时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}) =====")
        
        total_start_time = time.time()
        
        try:
            # 初始化流程
            self.initialize_pipeline()
            
            # 运行各个步骤
            self.run_probe_scoring()
            self.run_coverage_analysis()
            self.run_probe_mapping()
            
            total_elapsed_time = time.time() - total_start_time
            print(f"\n===== 玉米探针分析流程完成 (总耗时: {format_run_time(total_elapsed_time)}) =====")
            
            # 打印结果文件位置
            self.print_result_locations()
            
        except Exception as e:
            print(f"\n分析流程执行过程中出错: {str(e)}")
            import traceback
            traceback.print_exc()
            return False
        
        return True
    
    def print_result_locations(self):
        """打印结果文件位置"""
        print(f"\n===== 分析结果文件位置 ======")
        
        if 'score' in self.steps:
            print(f"探针评分结果: {self.config.OutputFiles.SCORE_RESULTS}")
            print(f"高特异性探针文件: {self.config.OutputFiles.HIGH_QUALITY_PROBES}")
        
        if 'coverage' in self.steps:
            print(f"覆盖分析报告: {self.config.OutputFiles.COVERAGE_REPORT}")
            print(f"覆盖分析图表: {self.config.OutputFiles.COVERAGE_PLOT}")
            print(f"分数分布图表: {self.config.OutputFiles.SCORE_DISTRIBUTION_PLOT}")
        
        if 'mapping' in self.steps:
            print(f"探针映射文本报告: {self.config.OutputFiles.MAPPING_RESULTS}")
            print(f"探针映射CSV报告: {self.config.OutputFiles.MAPPING_RESULTS_CSV}")
            print(f"探针映射HTML报告: {self.config.OutputFiles.MAPPING_RESULTS_HTML}")

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='玉米探针分析项目')
    
    # 添加步骤参数
    parser.add_argument('--steps', 
                      type=str, 
                      default='score,coverage,mapping', 
                      help='指定要运行的步骤，可选值为score,coverage,mapping，多个步骤用逗号分隔')
    
    # 添加配置文件参数
    parser.add_argument('--config', 
                      type=str, 
                      default=None, 
                      help='自定义配置文件路径')
    
    return parser.parse_args()

def main():
    """主函数"""
    # 解析命令行参数
    args = parse_arguments()
    
    # 解析步骤参数
    steps = [s.strip() for s in args.steps.split(',')]
    
    # 加载自定义配置（如果提供）
    custom_config = None
    if args.config:
        if os.path.exists(args.config):
            # 这里可以添加自定义配置加载逻辑
            print(f"警告: 自定义配置功能尚未实现，将使用默认配置")
        else:
            print(f"错误: 自定义配置文件不存在: {args.config}")
            sys.exit(1)
    
    # 创建并运行分析流程
    pipeline = ProbeAnalysisPipeline(steps=steps, custom_config=custom_config)
    success = pipeline.run()
    
    # 根据执行结果设置退出码
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()