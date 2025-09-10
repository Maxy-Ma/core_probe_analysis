#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通用工具函数模块
"""
import os
import sys
import logging
import chardet
from io import StringIO

# 尝试导入Biopython
try:
    from Bio import SeqIO
except ImportError:
    print("警告: 未找到Biopython库。请运行 'pip install biopython' 来安装。")

# 设置日志配置
def setup_logging(log_level='INFO'):
    """设置日志配置"""
    level_map = {
        'DEBUG': logging.DEBUG,
        'INFO': logging.INFO,
        'WARNING': logging.WARNING,
        'ERROR': logging.ERROR
    }
    
    level = level_map.get(log_level.upper(), logging.INFO)
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    return logging.getLogger(__name__)

# 自动检测文件编码并读取内容
def read_file_with_detect_encoding(file_path):
    """检测文件编码并读取内容"""
    try:
        # 先检测文件编码
        with open(file_path, 'rb') as f:
            raw_data = f.read(min(1024, os.path.getsize(file_path)))  # 只读取前1KB进行检测
            result = chardet.detect(raw_data)
            encoding = result['encoding'] or 'utf-8'  # 如果检测失败，默认使用utf-8
        
        # 尝试使用检测到的编码读取文件
        try:
            with open(file_path, 'r', encoding=encoding) as f:
                return f.read(), encoding
        except UnicodeDecodeError:
            # 如果使用检测到的编码失败，尝试使用utf-8并忽略错误
            print(f"警告: 使用检测到的编码 '{encoding}' 解码失败，尝试使用utf-8忽略错误")
            with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
                return f.read(), 'utf-8-replace'
    except Exception as e:
        print(f"读取文件时出错: {str(e)}")
        # 尝试创建一个清理后的版本
        try:
            clean_file = file_path.replace('.fasta', '_clean.fasta')
            with open(file_path, 'rb') as f_in, open(clean_file, 'wb') as f_out:
                for line in f_in:
                    clean_line = bytearray()
                    for b in line:
                        if b < 128 or b in (10, 13):  # 保留ASCII字符和换行符
                            clean_line.append(b)
                        else:
                            clean_line.append(63)  # '?'字符
                    f_out.write(clean_line)
            print(f"已创建清理后的文件: {clean_file}")
            
            # 从清理后的文件中读取
            with open(clean_file, 'r', encoding='ascii') as f:
                return f.read(), 'ascii-cleaned'
        except Exception as e2:
            raise Exception(f"创建清理后的文件并读取时也出错: {str(e2)}")

# 读取FASTA文件
def read_fasta_file(file_path):
    """读取FASTA文件并返回序列字典"""
    sequences = {}
    
    try:
        # 尝试使用BioPython读取
        content, encoding = read_file_with_detect_encoding(file_path)
        for record in SeqIO.parse(StringIO(content), "fasta"):
            seq_id = record.id
            seq = str(record.seq).upper()
            sequences[seq_id] = seq
        
        print(f"成功读取 {len(sequences)} 个序列 (编码: {encoding})")
        return sequences
    except Exception as e:
        print(f"使用BioPython读取FASTA文件时出错: {str(e)}")
        # 尝试手动解析FASTA格式
        print("尝试手动解析FASTA格式...")
        return parse_fasta_manually(file_path)

# 手动解析FASTA格式
def parse_fasta_manually(file_path):
    """手动解析FASTA格式文件"""
    sequences = {}
    
    try:
        content, encoding = read_file_with_detect_encoding(file_path)
        current_id = None
        current_seq = []
        
        for line in content.split('\n'):
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # 保存前一个序列
                if current_id:
                    sequences[current_id] = ''.join(current_seq).upper()
                    current_seq = []
                # 提取新的序列ID
                # 只取'>'后面的第一部分作为ID
                current_id = line[1:].split()[0] if ' ' in line[1:] else line[1:]
            else:
                # 累加序列数据
                if current_id:
                    current_seq.append(line)
        
        # 保存最后一个序列
        if current_id:
            sequences[current_id] = ''.join(current_seq).upper()
        
        print(f"手动解析成功读取 {len(sequences)} 个序列 (编码: {encoding})")
        return sequences
    except Exception as e:
        print(f"手动解析FASTA文件时出错: {str(e)}")
        return {}

# 写入FASTA文件
def write_fasta_file(sequences, output_path):
    """将序列字典写入FASTA文件"""
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n{seq}\n")
        print(f"成功写入 {len(sequences)} 个序列到 {output_path}")
        return True
    except Exception as e:
        print(f"写入FASTA文件时出错: {str(e)}")
        return False

# 读取BLAST输出文件
def read_blast_output(file_path):
    """读取BLAST输出文件并按探针ID组织结果"""
    blast_results = {}
    
    try:
        content, encoding = read_file_with_detect_encoding(file_path)
        for line in content.strip().split('\n'):
            if not line or line.startswith('#'):
                continue
            
            # 解析blast结果行
            parts = line.strip().split('\t')
            if len(parts) < 12:
                print(f"警告: blast行数据不足12列: {len(parts)}列")
                continue
            
            # 提取探针ID (第一列)
            probe_id = parts[0]
            
            # 将该行添加到对应探针的结果列表中
            if probe_id not in blast_results:
                blast_results[probe_id] = []
            blast_results[probe_id].append(line)
        
        print(f"成功读取BLAST结果，找到 {len(blast_results)} 个唯一探针 (编码: {encoding})")
        return blast_results
    except Exception as e:
        print(f"读取BLAST输出文件时出错: {str(e)}")
        return {}

# 检查文件是否存在
def check_file_exists(file_path):
    """检查文件是否存在"""
    if not os.path.exists(file_path):
        print(f"错误: 文件不存在: {file_path}")
        return False
    return True

# 确保目录存在
def ensure_dir_exists(dir_path):
    """确保目录存在，如果不存在则创建"""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(f"创建目录: {dir_path}")

# 格式化运行时间
def format_run_time(seconds):
    """将秒数格式化为人类可读的时间字符串"""
    if seconds < 60:
        return f"{seconds:.2f} 秒"
    elif seconds < 3600:
        minutes, seconds = divmod(seconds, 60)
        return f"{int(minutes)} 分 {int(seconds)} 秒"
    else:
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        return f"{int(hours)} 时 {int(minutes)} 分 {int(seconds)} 秒"

# 获取文件名（不含路径）
def get_filename_without_path(file_path):
    """获取文件名（不含路径）"""
    return os.path.basename(file_path)

# 获取文件名（不含扩展名）
def get_filename_without_extension(file_path):
    """获取文件名（不含扩展名）"""
    return os.path.splitext(os.path.basename(file_path))[0]

# 主函数（用于测试）
if __name__ == '__main__':
    # 测试文件读取功能
    if len(sys.argv) > 1:
        test_file = sys.argv[1]
        if os.path.exists(test_file):
            print(f"测试读取文件: {test_file}")
            content, encoding = read_file_with_detect_encoding(test_file)
            print(f"文件编码: {encoding}")
            print(f"文件内容前100个字符: {content[:100] if len(content) > 100 else content}")
        else:
            print(f"文件不存在: {test_file}")
    else:
        print("请提供测试文件路径")