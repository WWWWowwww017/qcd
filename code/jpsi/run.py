#!/usr/bin/env python3
"""可从任意工作目录运行，只依赖本包的父目录。"""
# 使导入不依赖调用者当前所在的工作目录。
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from jpsi.__main__ import main

if __name__ == "__main__":
    main()
