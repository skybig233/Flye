# -*- coding: utf-8 -*-
# @Time : 2021/4/28 13:27
# @Author : Jiangzhesheng
# @File : sam.py
# @Software: PyCharm
import sys
from typing import Iterator
from bioclass.file_object import File_object

QUERY_NAME=0
FLAG=1
REFERENCE_NAME=2
POSITION=3
MAP_QUALITY=4
CIGAR=5

class Sam:
    def __init__(self,line=''):
        list=line.split()
        try:
            self.query_name=list[QUERY_NAME]
            self.flag=list[FLAG]
            self.ref_name=list[REFERENCE_NAME]
            self.position=int(list[POSITION])
            self.mapq=int(list[MAP_QUALITY])
            self.cigar=list[CIGAR]
        except IndexError:raise TypeError(line,'is not a Sam')

    def __str__(self):
        tmp=[self.query_name,self.flag,self.ref_name,
             self.position,self.mapq,self.cigar]
        tmp=map(str,tmp)
        return '\t'.join(tmp)

class SamFile(File_object):

    def __init__(self, path: str):
        super().__init__(path)

    def __iter__(self,headerlines:int=0)->Iterator[Sam]:
        for i in super().__iter__(lines_per_unit=1,headerlines=headerlines):
            yield Sam(i)

class GeneralSamFile(SamFile):
    """
    一个通用的Sam文件类
    """
    def __init__(self, path: str)-> None :
        super().__init__(path)

    def __iter__(self)->Iterator[Sam]:
        """
        对SamFile迭代，产生Sam
        :return:一个装Sam的迭代器
        """
        # 首先判断sam文件是否有注释信息，以@开头，统计需要抛弃的注释行数
        headers=0
        with open(self.path) as f:
            for i in f:
                if i[:3] in ['@PG','@SQ','@HD','@RG','@CO']:
                    headers += 1
                else:
                    break

        # 再迭代
        # for i in super().__iter__(headerlines=headers):
        #     yield i
        return super().__iter__(headerlines=headers)

def main(argv):
    pass

if __name__ == '__main__':
    main(sys.argv)