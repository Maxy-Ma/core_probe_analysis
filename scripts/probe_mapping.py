#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
探针映射模块：将高特异性探针映射回PAV区域
"""
import os
import re
import time
from collections import defaultdict

# 导入工具函数
from utils import read_fasta_file, ensure_dir_exists

class ProbeMapper:
    """探针映射器类：用于将高特异性探针映射回PAV区域"""
    
    def __init__(self, config):
        """初始化探针映射器"""
        self.config = config
        self.high_quality_probes = {}
        self.pav_sequences = {}
        self.mapping_results = {}
    
    def load_high_quality_probes(self, probes_file=None):
        """加载高特异性探针"""
        file_path = probes_file or self.config.OutputFiles.HIGH_QUALITY_PROBES
        print(f"\n加载高特异性探针文件: {file_path}")
        
        # 检查文件是否存在
        if not os.path.exists(file_path):
            print(f"错误: 高特异性探针文件不存在: {file_path}")
            print("请先运行探针评分流程生成高特异性探针文件")
            return {}
        
        self.high_quality_probes = read_fasta_file(file_path)
        return self.high_quality_probes
    
    def load_pav_sequences(self, pav_file=None):
        """加载PAV序列"""
        file_path = pav_file or self.config.InputFiles.PAV_FILE
        print(f"\n加载PAV序列文件: {file_path}")
        self.pav_sequences = read_fasta_file(file_path)
        return self.pav_sequences
    
    def map_probes_to_pav(self):
        """将探针映射到PAV序列"""
        print(f"\n开始将探针映射到PAV序列...")
        start_time = time.time()
        
        # 初始化映射结果
        self.mapping_results = defaultdict(list)
        unmapped_probes = []
        
        # 将探针ID映射到pav序列ID
        for probe_id in self.high_quality_probes.keys():
            # 使用正则表达式提取基础序列ID
            # 探针ID格式如: 1_qLA2-1_upstream_1_pos0 -> 对应pav序列ID: 1_qLA2-1
            match = re.match(r'^([^_]+_[^_]+-\d+)_', probe_id)
            
            if match:
                pav_id = match.group(1)
                self.mapping_results[pav_id].append(probe_id)
            else:
                # 如果正则匹配失败，尝试简单的前缀匹配
                # 查找最接近的pav序列ID
                found = False
                for key in self.pav_sequences.keys():
                    if probe_id.startswith(key):
                        pav_id = key
                        self.mapping_results[pav_id].append(probe_id)
                        found = True
                        break
                
                if not found:
                    unmapped_probes.append(probe_id)
                    print(f"警告: 无法将探针 {probe_id} 映射到任何PAV序列ID")
        
        elapsed_time = time.time() - start_time
        
        # 统计映射结果
        total_probes = len(self.high_quality_probes)
        mapped_probes = sum(len(probes) for probes in self.mapping_results.values())
        mapped_pav_sequences = len(self.mapping_results)
        
        print(f"映射完成！")
        print(f"  处理时间: {elapsed_time:.2f} 秒")
        print(f"  总高特异性探针数: {total_probes}")
        print(f"  成功映射的探针数: {mapped_probes} ({mapped_probes/total_probes*100:.2f}%)")
        print(f"  未成功映射的探针数: {len(unmapped_probes)} ({len(unmapped_probes)/total_probes*100:.2f}%)")
        print(f"  找到对应的PAV序列数: {mapped_pav_sequences}")
        
        # 保存未映射的探针到文件
        if unmapped_probes:
            unmapped_file = os.path.join(self.config.RESULTS_DIR, 'unmapped_probes.txt')
            with open(unmapped_file, 'w', encoding='utf-8') as f:
                f.write("# 未成功映射到PAV序列的探针ID\n")
                for probe_id in unmapped_probes:
                    f.write(f"{probe_id}\n")
            print(f"  未映射的探针列表已保存到: {unmapped_file}")
        
        return self.mapping_results
    
    def generate_text_report(self, output_file=None):
        """生成文本格式的映射报告"""
        file_path = output_file or self.config.OutputFiles.MAPPING_RESULTS
        ensure_dir_exists(os.path.dirname(file_path))
        
        print(f"\n生成文本格式映射报告: {file_path}")
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write("# 高特异性探针与PAV序列的映射结果\n")
            f.write(f"# 高特异性分数阈值: {self.config.Params.SCORE_THRESHOLD}\n")
            f.write(f"# 共找到 {len(self.high_quality_probes)} 个高特异性探针，成功映射 {sum(len(probes) for probes in self.mapping_results.values())} 个\n")
            f.write(f"# 找到 {len(self.mapping_results)} 个对应的PAV序列\n")
            f.write("#\n")
            f.write("# 格式说明:\n")
            f.write("# 1. 每个PAV序列ID后面跟着关联的探针列表\n")
            f.write("# 2. 探针按名称排序\n")
            f.write("# 3. 显示PAV序列的长度和前50个碱基\n")
            f.write("#\n\n")
            
            # 按PAV序列ID排序输出
            for pav_id in sorted(self.mapping_results.keys()):
                probes = self.mapping_results[pav_id]
                # 对探针进行排序
                probes.sort()
                
                # 获取PAV序列信息
                seq = self.pav_sequences.get(pav_id, "")
                seq_len = len(seq)
                seq_preview = seq[:50] + "..." if len(seq) > 50 else seq
                
                # 输出PAV序列信息和关联的探针
                f.write(f"====================================================\n")
                f.write(f"PAV序列ID: {pav_id}\n")
                f.write(f"序列长度: {seq_len} bp\n")
                f.write(f"序列前50bp: {seq_preview}\n")
                f.write(f"关联的高特异性探针数量: {len(probes)}\n")
                f.write(f"----------------------------------------------------\n")
                
                # 输出每个关联的探针
                for i, probe_id in enumerate(probes):
                    f.write(f"探针 {i+1}: {probe_id}\n")
                
                f.write("====================================================\n\n")
        
        print(f"文本格式映射报告已保存到 {file_path}")
        return file_path
    
    def generate_csv_report(self, output_file=None):
        """生成CSV格式的映射报告"""
        file_path = output_file or self.config.OutputFiles.MAPPING_RESULTS_CSV
        ensure_dir_exists(os.path.dirname(file_path))
        
        print(f"\n生成CSV格式映射报告: {file_path}")
        
        with open(file_path, 'w', encoding='utf-8', newline='') as f:
            # 写入CSV标题行
            f.write("PAV序列ID,探针ID,探针数量,序列长度\n")
            
            # 按PAV序列ID排序输出
            for pav_id in sorted(self.mapping_results.keys()):
                probes = self.mapping_results[pav_id]
                seq = self.pav_sequences.get(pav_id, "")
                seq_len = len(seq)
                
                # 为每个探针写入一行
                for probe_id in probes:
                    f.write(f"{pav_id},{probe_id},{len(probes)},{seq_len}\n")
        
        print(f"CSV格式映射报告已保存到 {file_path}")
        return file_path
    
    def generate_html_report(self, output_file=None):
        """生成HTML格式的映射报告"""
        file_path = output_file or self.config.OutputFiles.MAPPING_RESULTS_HTML
        ensure_dir_exists(os.path.dirname(file_path))
        
        print(f"\n生成HTML格式映射报告: {file_path}")
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write('<!DOCTYPE html>\n')
            f.write('<html lang="zh-CN">\n')
            f.write('<head>\n')
            f.write('    <meta charset="UTF-8">\n')
            f.write('    <meta name="viewport" content="width=device-width, initial-scale=1.0">\n')
            f.write('    <title>高特异性探针与PAV序列映射报告</title>\n')
            f.write('    <style>\n')
            f.write('        body { font-family: Arial, sans-serif; margin: 20px; }\n')
            f.write('        h1 { color: #333; }\n')
            f.write('        .summary { background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px; }\n')
            f.write('        .sequence-container { border: 1px solid #ddd; border-radius: 5px; margin-bottom: 20px; overflow: hidden; }\n')
            f.write('        .sequence-header { background-color: #4CAF50; color: white; padding: 10px; }\n')
            f.write('        .sequence-body { padding: 15px; }\n')
            f.write('        .probe-list { list-style-type: none; padding: 0; }\n')
            f.write('        .probe-list li { padding: 5px 0; border-bottom: 1px solid #eee; }\n')
            f.write('        .probe-list li:last-child { border-bottom: none; }\n')
            f.write('        .footer { margin-top: 30px; text-align: center; color: #666; font-size: 0.9em; }\n')
            f.write('        .search-container { margin-bottom: 20px; }\n')
            f.write('        #search-input { padding: 8px; width: 300px; }\n')
            f.write('        .pagination { margin-top: 20px; text-align: center; }\n')
            f.write('        .pagination a { margin: 0 5px; text-decoration: none; color: #4CAF50; }\n')
            f.write('        .pagination a.active { font-weight: bold; }\n')
            f.write('    </style>\n')
            f.write('</head>\n')
            f.write('<body>\n')
            f.write('    <h1>高特异性探针与PAV序列映射报告</h1>\n')
            f.write('    <div class="summary">\n')
            f.write(f'        <p>高特异性分数阈值: {self.config.Params.SCORE_THRESHOLD}</p>\n')
            f.write(f'        <p>找到的高特异性探针数量: {len(self.high_quality_probes)}</p>\n')
            f.write(f'        <p>成功映射的探针数量: {sum(len(probes) for probes in self.mapping_results.values())}</p>\n')
            f.write(f'        <p>对应的PAV序列数量: {len(self.mapping_results)}</p>\n')
            f.write('    </div>\n')
            f.write('    <div class="search-container">\n')
            f.write('        <input type="text" id="search-input" placeholder="搜索PAV序列ID或探针ID...">\n')
            f.write('    </div>\n')
            f.write('    <h2>映射结果详情</h2>\n')
            f.write('    <div id="results-container">\n')
            
            # 按PAV序列ID排序输出
            for pav_id in sorted(self.mapping_results.keys()):
                probes = self.mapping_results[pav_id]
                probes.sort()
                seq = self.pav_sequences.get(pav_id, "")
                seq_len = len(seq)
                seq_preview = seq[:50] + "..." if len(seq) > 50 else seq
                
                f.write('    <div class="sequence-container">\n')
                f.write(f'        <div class="sequence-header">\n')
                f.write(f'            <h3>PAV序列ID: {pav_id}</h3>\n')
                f.write('        </div>\n')
                f.write('        <div class="sequence-body">\n')
                f.write(f'            <p><strong>序列长度:</strong> {seq_len} bp</p>\n')
                f.write(f'            <p><strong>序列前50bp:</strong> {seq_preview}</p>\n')
                f.write(f'            <p><strong>关联的高特异性探针数量:</strong> {len(probes)}</p>\n')
                f.write('            <ul class="probe-list">\n')
                for probe_id in probes:
                    f.write(f'                <li>{probe_id}</li>\n')
                f.write('            </ul>\n')
                f.write('        </div>\n')
                f.write('    </div>\n')
            
            f.write('    </div>\n')
            f.write('    <div class="footer">\n')
            f.write(f'        <p>报告生成时间: {time.strftime("%Y-%m-%d %H:%M:%S")}</p>\n')
            f.write('    </div>\n')
            f.write('    <script>\n')
            f.write('        // 搜索功能\n')
            f.write('        document.getElementById("search-input").addEventListener("input", function() {\n')
            f.write('            const searchTerm = this.value.toLowerCase();\n')
            f.write('            const containers = document.querySelectorAll(".sequence-container");\n')
            f.write('            let visibleCount = 0;\n')
            f.write('            \n')
            f.write('            containers.forEach(container => {\n')
            f.write('                const text = container.textContent.toLowerCase();\n')
            f.write('                if (text.includes(searchTerm)) {\n')
            f.write('                    container.style.display = "block";\n')
            f.write('                    visibleCount++;\n')
            f.write('                } else {\n')
            f.write('                    container.style.display = "none";\n')
            f.write('                }\n')
            f.write('            });\n')
            f.write('            \n')
            f.write('            // 如果没有匹配结果，显示提示\n')
            f.write('            const resultsContainer = document.getElementById("results-container");\n')
            f.write('            let noResultsMsg = document.getElementById("no-results-msg");\n')
            f.write('            \n')
            f.write('            if (visibleCount === 0) {\n')
            f.write('                if (!noResultsMsg) {\n')
            f.write('                    noResultsMsg = document.createElement("div");\n')
            f.write('                    noResultsMsg.id = "no-results-msg";\n')
            f.write('                    noResultsMsg.textContent = "未找到匹配的结果。";\n')
            f.write('                    resultsContainer.appendChild(noResultsMsg);\n')
            f.write('                }\n')
            f.write('                noResultsMsg.style.display = "block";\n')
            f.write('            } else if (noResultsMsg) {\n')
            f.write('                noResultsMsg.style.display = "none";\n')
            f.write('            }\n')
            f.write('        });\n')
            f.write('    </script>\n')
            f.write('</body>\n')
            f.write('</html>\n')
        
        print(f"HTML格式映射报告已保存到 {file_path}")
        return file_path
    
    def run(self):
        """运行完整的探针映射流程"""
        print("===== 探针映射流程开始 =====")
        
        # 加载数据
        self.load_high_quality_probes()
        self.load_pav_sequences()
        
        # 如果没有高特异性探针，跳过映射步骤
        if not self.high_quality_probes:
            print("警告: 没有高特异性探针可供映射")
            return None
        
        # 执行映射
        self.map_probes_to_pav()
        
        # 生成报告
        self.generate_text_report()
        self.generate_csv_report()
        self.generate_html_report()
        
        print("===== 探针映射流程完成 =====")
        
        return self.mapping_results

# 主函数（用于测试）
if __name__ == '__main__':
    # 导入配置
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from config import config, create_directories
    
    # 创建必要的目录
    create_directories()
    
    # 创建映射器实例并运行
    mapper = ProbeMapper(config)
    mapper.run()