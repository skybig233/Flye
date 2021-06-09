# -*- coding: utf-8 -*-
# @Time : 2020/9/27 9:42
# @Author : Jiangzhesheng
# @File : file_object.py
# @Software: PyCharm
import sys
import os
from typing import Iterator
class File_object:
    # """
    # 在生信的某些文件中，一个信息单元（unit）的信息，由文件多行提供
    # 例如fastq文件里面的fastq，一个信息单元（fastq），由4行信息构成：id、base、orient、quality
    # File_object是这些文件的父类
    # 属性：path，文件的路径
    # 方法：iter，迭代器，每次调用返回一个信息单元（unit）
    # """
    def __init__(self,path:str):
        self.path=path
        super().__init__()

    def __iter__(self, lines_per_unit:int,headerlines:int)->Iterator[str]:
        """
        :param lines_per_unit: 一个生信单元用多少行描述
        :param headerlines:文件头部，不用于描述生信单元的行数
        :return: 文件中的生信单元
        """

        with open(self.path,mode='r') as file:
            # 先舍弃无用行
            for i in range(headerlines):
                file.readline()
            while True:
                # 再读入每个单元
                unit=''
                for i in range(lines_per_unit):
                    unit+=file.readline()
                if not unit:
                    break
                yield unit

def main(argv):
    pass

if __name__ == '__main__':
    main(sys.argv)
