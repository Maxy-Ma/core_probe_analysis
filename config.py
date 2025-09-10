#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置文件：导入配置类并创建全局配置对象
"""
# 导入配置类
from config_class import config, create_directories, validate_input_files

# 这个文件只是一个简单的导入文件，实际的配置定义在config_class.py中
# 这样设计的目的是将配置的定义和使用分离，便于维护和扩展

# 可以在这里添加任何项目级别的全局变量或常量
# 例如：
# PROJECT_VERSION = '1.0.0'
# PROJECT_NAME = '玉米探针分析项目'

# 导出主要的配置对象和函数
__all__ = ['config', 'create_directories', 'validate_input_files']