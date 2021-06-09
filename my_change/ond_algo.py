# -*- coding: utf-8 -*-
# @Time : 2020/11/24 10:47
# @Author : Jiangzhesheng
# @File : algo.py
# @Software: PyCharm
import sys
import time
class Point:
    def __init__(self,x,k) -> None:
        self.x=x
        self.k=k
        self.y=x-k

    def __str__(self) -> str:
        return '(%s,%s)' % (self.x,self.y)

def list_get(index:int, list:[], default=None):
    return list[index] if 0<= index < len(list) else default

def ond_algo(A:str, B:str)->[[]]:
    """
    改进后的OND算法，IDS权值都为1，并且返回一个v数组，后续根据v数组计算D值和路径（编辑脚本）
    :param A:
    :param B:
    :return:v数组
    """
    n=len(A)
    m=len(B)
    max_value=m+n
    v=[]
    x,y=0,0
    while x < n and y < m and A[x] == B[y]:
        x, y = x + 1, y + 1
    v.append([x])
    # 边界值处理。。。。例如('aaaa','aaaa')、('','')
    if x>=n and y>=m:
        return v
    for d in range(1,max_value):
        tmp_v=[-1]*(2*d+1)
        for k in range(-d,d+1,1):
            x1=list_get(list=v[d-1],index=k+d-2,default=-1)
            x2=list_get(list=v[d-1],index=k+d-1,default=-1)
            x3=list_get(list=v[d-1],index=k+d,default=-1)
            p1 = Point(x=x1, k=k-1)#横着走的点
            p2= Point(x=x2,k=k)#substitution
            p3 = Point(x=x3, k=k+1)#竖着走的点
            x=max([p1.x+1,p2.x+1,p3.x])
            y=x-k
            while x<n and y<m and A[x]==B[y]:
                x,y=x+1,y+1
            tmp_v[k+d]=x
            if x>=n and y>=m:
                v.append(tmp_v)
                return v
        v.append(tmp_v)

def get_edit_distance(v:[[]])->int:
    return len(v)-1

def get_edit_script(v:[[]])->[]:
    script=[]
    d=len(v)-1
    k=v[-1].index(-1)-d-1 if -1 in v[-1] else d
    for d in range(len(v)-1,0,-1):
        x1 = list_get(list=v[d - 1], index=k + d - 2, default=-1)
        x2 = list_get(list=v[d - 1], index=k + d - 1, default=-1)
        x3 = list_get(list=v[d - 1], index=k + d, default=-1)
        p1 = Point(x=x1, k=k - 1)  # 横着走的点
        p2 = Point(x=x2, k=k)  # substitution
        p3 = Point(x=x3, k=k + 1)  # 竖着走的点
        p_list=[(p1.x+1,p1),(p2.x+1,p2),(p3.x,p3)]
        p_list.sort(key=lambda p:p[0])
        if p_list[-1][1]==p1:
            script.insert(0,(p_list[-1][1].x,'D'))
        elif p_list[-1][1]==p2:
            script.insert(0,(p_list[-1][1].x,'S'+str(p_list[-1][1].y)))
        elif p_list[-1][1]==p3:
            script.insert(0,(p_list[-1][1].x,'I'+str(p_list[-1][1].y)))
        k=p_list[-1][1].k
    return script

def main(argv):
    pass

if __name__ == '__main__':
    main(sys.argv)