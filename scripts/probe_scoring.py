#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
探针评分模块：实现使用改进算法计算探针特异性分数
"""
import os
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

# 导入工具函数
from utils import read_fasta_file, read_blast_output, ensure_dir_exists, format_run_time

class ProbeScorer:
    """探针评分器类：用于计算探针特异性分数"""
    
    def __init__(self, config):
        """初始化探针评分器"""
        self.config = config
        self.probes = {}
        self.blast_results = {}
        self.scores = {}
        
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
        
        self.blast_results = defaultdict(list)
        
        for file_path in files_to_load:
            print(f"  正在加载: {os.path.basename(file_path)}")
            file_results = read_blast_output(file_path)
            
            # 将当前文件的结果合并到总结果中
            for probe_id, hits in file_results.items():
                self.blast_results[probe_id].extend(hits)
        
        print(f"成功加载 {len(self.blast_results)} 个唯一探针的BLAST结果")
        return self.blast_results
    
    def calculate_probe_score(self, probe_id, hits):
        """
        计算单个探针的特异性分数（改进算法）
        
        步骤3: 计算原始非特异性指数 (Raw Non-Specificity Index)
        Raw_Non_Specificity = Σ_{i=1}^{N} (nident_i / L)
        
        步骤4: 计算原始特异性分数 (Raw Specificity Score)
        Raw_Specificity_Score = 1 / (1 + Raw_Non_Specificity)
        
        步骤5: 标准化与最终分数 (Normalization & Final Score)
        Final_Specificity_Score = 100 * Raw_Specificity_Score
        """
        try:
            # 获取探针序列，如果不存在则使用空字符串
            probe_seq = self.probes.get(probe_id, "")
            probe_length = len(probe_seq)
            
            # 初始化参数
            raw_non_specificity = 0.0
            total_hits = 0
            perfect_matches = 0
            
            # 收集染色体信息
            matched_chromosomes = set()
            
            for hit in hits:
                # 解析blast结果行
                parts = hit.strip().split('\t')
                if len(parts) < 12:
                    continue
                
                # 获取染色体信息（第2列）
                chromosome = parts[1]
                matched_chromosomes.add(chromosome)
                
                # 获取相似度和比对长度信息
                try:
                    # 步骤3: 计算原始非特异性指数
                    # Raw_Non_Specificity = Σ_{i=1}^{N} (nident_i / L)
                    pident = float(parts[2])  # 相似度百分比
                    length = int(parts[3])  # 比对长度
                    
                    # 计算匹配碱基数
                    nident_i = (pident / 100.0) * length
                    
                    # 累加非特异性指数
                    if probe_length > 0:
                        raw_non_specificity += (nident_i / probe_length)
                    
                    # 统计完全匹配
                    if pident == 100.0 and int(parts[4]) == 0:  # 100%相似度且不错配
                        perfect_matches += 1
                    
                    total_hits += 1
                except (ValueError, IndexError) as e:
                    print(f"解析blast数据时出错: {str(e)}")
                    continue
            
            # 步骤4: 计算原始特异性分数
            raw_specificity_score = 1.0 / (1.0 + raw_non_specificity)
            
            # 步骤5: 标准化与最终分数
            final_score = 100.0 * raw_specificity_score
            
            # 特殊情况处理
            if total_hits == 0:
                final_score = 100.0  # 没有找到任何匹配，给满分100分
            
            if not probe_seq:
                final_score = min(final_score, 30.0)  # 探针序列为空，最多给30分
                print(f"  警告: 探针 {probe_id} 在FASTA文件中未找到，分数上限为30分")
            
            return {
                'probe_id': probe_id,
                'score': round(final_score, 2),
                'raw_non_specificity': round(raw_non_specificity, 4),
                'perfect_matches': perfect_matches,
                'total_hits': total_hits,
                'matched_chromosomes': ','.join(sorted(matched_chromosomes)),
                'probe_length': probe_length,
                'in_fasta': probe_id in self.probes
            }
        except Exception as e:
            print(f"计算探针 {probe_id} 分数时出错: {str(e)}")
            return None
    
    def calculate_probes_batch(self, probe_ids_hits):
        """批量计算探针分数（用于多进程处理）"""
        results = []
        for probe_id, hits in probe_ids_hits:
            result = self.calculate_probe_score(probe_id, hits)
            if result:
                results.append(result)
        return results

    @staticmethod
    def static_calculate_probes_batch(probes, probe_ids_hits):
        """静态方法版本的批量计算探针分数（专门用于多进程处理）"""
        results = []
        for probe_id, hits in probe_ids_hits:
            # 手动模拟calculate_probe_score方法的核心逻辑
            try:
                # 获取探针序列
                probe_seq = probes.get(probe_id, "")
                probe_length = len(probe_seq)
                
                # 初始化参数
                raw_non_specificity = 0.0
                total_hits = 0
                perfect_matches = 0
                matched_chromosomes = set()
                
                for hit in hits:
                    parts = hit.strip().split('\t')
                    if len(parts) < 12:
                        continue
                    
                    # 获取染色体信息
                    try:
                        chromosome = parts[1]
                        matched_chromosomes.add(chromosome)
                        
                        # 计算非特异性指数
                        pident = float(parts[2])
                        length = int(parts[3])
                        nident_i = (pident / 100.0) * length
                        
                        if probe_length > 0:
                            raw_non_specificity += (nident_i / probe_length)
                        
                        # 统计完全匹配
                        if pident == 100.0 and int(parts[4]) == 0:
                            perfect_matches += 1
                        
                        total_hits += 1
                    except (ValueError, IndexError):
                        continue
                
                # 计算分数
                raw_specificity_score = 1.0 / (1.0 + raw_non_specificity)
                final_score = 100.0 * raw_specificity_score
                
                # 特殊情况处理
                if total_hits == 0:
                    final_score = 100.0  # 没有找到任何匹配，给满分100分
                
                if not probe_seq:
                    final_score = min(final_score, 30.0)
                
                results.append({
                    'probe_id': probe_id,
                    'score': round(final_score, 2),
                    'raw_non_specificity': round(raw_non_specificity, 4),
                    'perfect_matches': perfect_matches,
                    'total_hits': total_hits,
                    'matched_chromosomes': ','.join(sorted(matched_chromosomes)) if matched_chromosomes else '',
                    'probe_length': probe_length,
                    'in_fasta': probe_id in probes
                })
            except Exception as e:
                print(f"计算探针 {probe_id} 分数时出错: {str(e)}")
                continue
        
        return results
    
    def calculate_all_scores(self):
        """计算所有探针的分数"""
        print(f"\n开始计算探针特异性分数...")
        start_time = time.time()
        
        # 准备探针列表
        probe_list = list(self.blast_results.items())
        total_probes = len(probe_list)
        print(f"  共有 {total_probes} 个探针需要计算分数")
        
        # 将探针分成批次以减少进程创建开销
        chunk_size = self.config.Params.CHUNK_SIZE
        batches = [probe_list[i:i+chunk_size] for i in range(0, total_probes, chunk_size)]
        
        # 使用多进程并行处理
        max_workers = self.config.Params.MAX_WORKERS
        print(f"  使用多进程处理 (进程数: {max_workers}, 每批: {chunk_size})")
        
        self.scores = []
        processed_count = 0
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务 - 使用静态方法避免序列化问题
            # 将探针字典和批次数据一起传入
            future_to_batch = {executor.submit(ProbeScorer.static_calculate_probes_batch, self.probes, batch): i for i, batch in enumerate(batches)}
            
            # 收集结果
            for future in as_completed(future_to_batch):
                batch_idx = future_to_batch[future]
                try:
                    batch_results = future.result()
                    self.scores.extend(batch_results)
                    processed_count += len(batch_results)
                    
                    # 显示进度
                    progress = min(int(processed_count / total_probes * 100), 100)
                    print(f"  处理进度: {processed_count}/{total_probes} ({progress}%)")
                except Exception as e:
                    print(f"  处理批次 {batch_idx} 时出错: {str(e)}")
        
        # 按分数排序
        self.scores.sort(key=lambda x: x['score'], reverse=True)
        
        elapsed_time = time.time() - start_time
        print(f"分数计算完成！")
        print(f"  处理时间: {format_run_time(elapsed_time)}")
        print(f"  成功计算 {len(self.scores)} 个探针的分数")
        
        # 生成统计信息
        self._generate_score_statistics()
        
        return self.scores
    
    def _generate_score_statistics(self):
        """生成分数统计信息"""
        if not self.scores:
            return
        
        # 计算平均分数
        avg_score = sum(r['score'] for r in self.scores) / len(self.scores)
        
        # 统计不同分数区间的探针数量
        high_score_probes = sum(1 for r in self.scores if r['score'] >= 80)
        medium_score_probes = sum(1 for r in self.scores if 60 <= r['score'] < 80)
        low_score_probes = sum(1 for r in self.scores if r['score'] < 60)
        
        # 统计在FASTA中找到的探针数量
        found_in_fasta = sum(1 for r in self.scores if r['in_fasta'])
        
        print(f"\n===== 分数统计信息 =====")
        print(f"平均分数: {round(avg_score, 2)}")
        print(f"高特异性探针(>=80分): {high_score_probes} ({high_score_probes/len(self.scores)*100:.2f}%)")
        print(f"中等特异性探针(60-79分): {medium_score_probes} ({medium_score_probes/len(self.scores)*100:.2f}%)")
        print(f"低特异性探针(<60分): {low_score_probes} ({low_score_probes/len(self.scores)*100:.2f}%)")
        print(f"在FASTA中找到的探针: {found_in_fasta} ({found_in_fasta/len(self.scores)*100:.2f}%)")
        
    def save_scores_to_file(self, output_file=None):
        """保存分数结果到文件"""
        file_path = output_file or self.config.OutputFiles.SCORE_RESULTS
        ensure_dir_exists(os.path.dirname(file_path))
        
        print(f"\n保存分数结果到文件: {file_path}")
        
        with open(file_path, 'w', encoding='utf-8') as f:
            # 写入标题行
            f.write("探针ID\t特异性分数\t原始非特异性指数\t完全匹配数\t总匹配数\t匹配染色体\t探针长度\t在FASTA中找到\n")
            
            # 写入每个探针的分数
            for result in self.scores:
                found_in_fasta = "是" if result['in_fasta'] else "否"
                f.write(
                    f"{result['probe_id']}\t"
                    f"{result['score']}\t"
                    f"{result['raw_non_specificity']}\t"
                    f"{result['perfect_matches']}\t"
                    f"{result['total_hits']}\t"
                    f"{result['matched_chromosomes']}\t"
                    f"{result['probe_length']}\t"
                    f"{found_in_fasta}\n"
                )
        
        print(f"成功保存 {len(self.scores)} 个探针的分数结果")
        return file_path
    
    def extract_high_quality_probes(self, output_file=None, threshold=None):
        """提取高特异性探针"""
        # 使用60分作为默认阈值
        score_threshold = threshold or 60.0
        file_path = output_file or self.config.OutputFiles.HIGH_QUALITY_PROBES
        
        print(f"\n提取高特异性探针（阈值: {score_threshold}分）")
        
        # 筛选高特异性探针
        high_quality_probes = {}
        for result in self.scores:
            if result['score'] >= score_threshold and result['in_fasta']:
                probe_id = result['probe_id']
                high_quality_probes[probe_id] = self.probes[probe_id]
        
        # 检查是否没有足够质量的探针
        is_low_quality_fallback = False
        if not high_quality_probes and self.scores and any(r['in_fasta'] for r in self.scores):
            # 没有探针分数高于阈值，取最高分的探针
            best_probe = None
            for result in self.scores:
                if result['in_fasta'] and (best_probe is None or result['score'] > best_probe['score']):
                    best_probe = result
            
            if best_probe:
                # 标记为质量不好的探针
                probe_id = best_probe['probe_id']
                high_quality_probes[probe_id] = self.probes[probe_id]
                is_low_quality_fallback = True
                print(f"  警告: 没有探针分数高于 {score_threshold} 分，取最高分探针作为备选")
                print(f"  最高分探针: {probe_id}，分数: {best_probe['score']} 分")
        
        # 保存到文件
        ensure_dir_exists(os.path.dirname(file_path))
        
        with open(file_path, 'w', encoding='utf-8') as f:
            # 写入文件头信息
            if is_low_quality_fallback:
                f.write(f"> # WARNING: LOW QUALITY PROBE (NO PROBES ABOVE {score_threshold} SCORE)\n")
                f.write(f"> # This is the highest scoring probe but may have poor specificity\n")
            else:
                f.write(f"> # High Quality Probes (Score >= {score_threshold})\n")
            f.write(f"> # Total: {len(high_quality_probes)} probes\n\n")
            
            # 写入探针序列，为低质量备选探针添加特殊标记
            for probe_id, seq in high_quality_probes.items():
                if is_low_quality_fallback:
                    f.write(f">{probe_id} [LOW_QUALITY_FALLBACK]\n{seq}\n")
                else:
                    f.write(f">{probe_id}\n{seq}\n")
        
        print(f"成功提取 {len(high_quality_probes)} 个高特异性探针到 {file_path}")
        return high_quality_probes, file_path
    
    def run(self):
        """运行完整的评分流程"""
        print("===== 探针评分流程开始 =====")
        
        # 加载数据
        self.load_probes()
        self.load_blast_results()
        
        # 计算分数
        self.calculate_all_scores()
        
        # 保存结果
        self.save_scores_to_file()
        
        # 提取高特异性探针
        self.extract_high_quality_probes()
        
        print("===== 探针评分流程完成 =====")
        
        return self.scores

# 主函数（用于测试）
if __name__ == '__main__':
    # 导入配置
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from config import config, create_directories
    
    # 创建必要的目录
    create_directories()
    
    # 创建评分器实例并运行
    scorer = ProbeScorer(config)
    scorer.run()