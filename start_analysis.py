#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
玉米探针分析项目 - Python启动脚本

使用方法：
1. 确保所有必要的输入文件都已放入data目录
2. 双击运行此脚本，或在命令提示符中运行：python start_analysis.py
"""
import os
import sys
import time
from datetime import datetime

# 添加项目根目录到Python路径
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def print_header():
    """打印程序头部信息"""
    print("=" * 80)
    print("                    玉米探针分析项目")
    print("=" * 80)
    print(f"启动时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"当前Python版本: {sys.version.split()[0]}")
    print(f"项目路径: {os.path.dirname(os.path.abspath(__file__))}")
    print("=" * 80)
    print()

def check_python_environment():
    """检查Python环境"""
    print("[步骤1/3] 检查Python环境...")
    
    # 检查Python版本
    major, minor = sys.version_info[:2]
    if major < 3 or (major == 3 and minor < 6):
        print("错误: 需要Python 3.6或更高版本！")
        print(f"当前版本: {major}.{minor}")
        return False
    
    print(f"✓ Python版本检查通过: {major}.{minor}")
    
    # 检查必要的Python库
    required_libs = ['Bio', 'matplotlib', 'pandas', 'numpy', 'chardet']
    missing_libs = []
    
    for lib in required_libs:
        try:
            __import__(lib)
            print(f"✓ 已安装库: {lib}")
        except ImportError:
            missing_libs.append(lib)
            print(f"✗ 未安装库: {lib}")
    
    if missing_libs:
        print("\n错误: 缺少必要的Python库！")
        print("请运行以下命令安装缺失的库：")
        print(f"pip install {' '.join(missing_libs)}")
        return False
    
    print("[步骤1/3] Python环境检查通过！\n")
    return True

def check_project_structure():
    """检查项目目录结构"""
    print("[步骤2/3] 检查项目目录结构...")
    
    # 获取项目根目录
    project_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 检查必要的子目录
    required_dirs = ['data', 'scripts', 'results', 'reports']
    missing_dirs = []
    
    for dir_name in required_dirs:
        dir_path = os.path.join(project_dir, dir_name)
        if not os.path.exists(dir_path):
            missing_dirs.append(dir_path)
            try:
                os.makedirs(dir_path)
                print(f"✓ 创建缺失的目录: {dir_path}")
            except Exception as e:
                print(f"✗ 创建目录失败: {dir_path}, 错误: {str(e)}")
        else:
            print(f"✓ 目录存在: {dir_path}")
    
    # 检查scripts目录下的必要文件
    scripts_dir = os.path.join(project_dir, 'scripts')
    required_scripts = ['main.py', 'config.py', 'utils.py', 'probe_scoring.py', 'coverage_analysis.py', 'probe_mapping.py']
    missing_scripts = []
    
    for script in required_scripts:
        script_path = os.path.join(scripts_dir, script)
        if not os.path.exists(script_path):
            missing_scripts.append(script_path)
            print(f"✗ 文件不存在: {script_path}")
        else:
            print(f"✓ 文件存在: {script}")
    
    if missing_scripts:
        print("\n错误: 缺少必要的脚本文件！")
        return False
    
    print("[步骤2/3] 项目目录结构检查通过！\n")
    return True

def run_analysis():
    """运行分析程序"""
    print("[步骤3/3] 开始运行分析程序...")
    
    # 添加scripts目录到Python路径
    project_dir = os.path.dirname(os.path.abspath(__file__))
    scripts_dir = os.path.join(project_dir, 'scripts')
    sys.path.append(scripts_dir)
    
    try:
        # 导入主模块
        import main
        
        print("✓ 成功导入主模块")
        print("\n开始执行分析流程，请稍候...\n")
        
        # 记录开始时间
        start_time = time.time()
        
        # 调用主函数
        main.main()
        
        # 记录结束时间
        elapsed_time = time.time() - start_time
        
        print(f"\n✓ 分析流程执行完成！")
        print(f"总耗时: {time.strftime('%H小时%M分钟%S秒', time.gmtime(elapsed_time))}")
        print(f"结果文件保存在: {os.path.join(project_dir, 'results')} 和 {os.path.join(project_dir, 'reports')}")
        
        return True
        
    except ImportError as e:
        print(f"✗ 导入模块失败: {str(e)}")
        print("请检查项目文件是否完整。")
        return False
    except Exception as e:
        print(f"✗ 分析过程中出现错误: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """主函数"""
    # 打印头部信息
    print_header()
    
    # 检查Python环境
    if not check_python_environment():
        print("\n环境检查失败，无法继续执行分析。")
        input("\n按回车键退出...")
        sys.exit(1)
    
    # 检查项目结构
    if not check_project_structure():
        print("\n项目结构检查失败，无法继续执行分析。")
        input("\n按回车键退出...")
        sys.exit(1)
    
    # 运行分析
    success = run_analysis()
    
    # 等待用户按键退出
    if success:
        input("\n分析已完成。按回车键退出...")
    else:
        input("\n分析过程中出现错误。按回车键退出...")

if __name__ == '__main__':
    main()