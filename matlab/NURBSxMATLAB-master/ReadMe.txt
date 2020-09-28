A2.1
A2.2
A2.3
A3.5 B-spline曲面 done A3.6 偏导矢
A4.3 NURBS曲面    done A4.4 偏导矢

A2.2 u= 4:VectorValue=0 (推导当knotVector中有重复节点值(在中间)出现时的公式 special case?)


1.rhino中knot数量少2. (definition: knots = CvPts + degree +1)
  解决：在knot向量最前面和最后面各加上一个knot，值同相邻knot。
2.当u属于[ui,ui+1)，u=ui+1时如何处理？(书上没看到)


IGES格式 / STEP格式
ShrinkTrimmedSrfToEdge 缩回已修剪曲面

p99 4.20 / S计算顺序 / 看p101 算法4.4 / 尝试表驱动
边界处考虑