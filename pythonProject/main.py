# # # 这是一个示例 Python 脚本。
# #
# # # 按 Shift+F10 执行或将其替换为您的代码。
# # # 按 双击 Shift 在所有地方搜索类、文件、工具窗口、操作和设置。
# #
# #
# # def print_hi(name):
# #     # 在下面的代码行中使用断点来调试脚本。
# #     print(f'Hi, {name}')  # 按 Ctrl+F8 切换断点。
# #
# #
# # # 按间距中的绿色按钮以运行脚本。
# # if __name__ == '__main__':
# #     # name = 'notify'
# #     # print_hi('PyCharm')
# #     # a = input() ctrl+/快速注释
# #     # a = int(input())
# #     # print(a)
# #     # print(type(a))
# #     stu = 6191611025
# #     name = 'tt'
# #     print('stu number is {} and name is {}'.format(stu, name))
# #     print(f'stu number is {stu} and name is {name}')
# # # 访问 https://www.jetbrains.com/help/pycharm/ 获取 PyCharm 帮助
# """
# import math
#
# print(math.ceil(2.1)) #向上取整
# print(math.floor(2.8)) #向下取整
# print(math.fabs(-2.5)) #绝对值
# print(math.pow(2, 3)) #2的3次方
# print(math.sqrt(4)) #开平方
#
# """
# import random
#
# print(random.random())  # [0,1)随机数
# print(random.randint(3, 6))  # [3,6]随机整数
# print(random.randrange(1, 10, 3))  # [1,10)区间内步长为3的随机整数
# print(random.uniform(1, 10))  # [1,10]随机浮点数
# a = range(10)
# print(random.choice(a))  # 从列表随机选取元素 range(10)=[0,10)整数
# b = range(1, 10)
# c=list(b)
# print(c)
# print(type(b))
# print(type(c))
# print(random.shuffle(c))  #直接打乱原有列表顺序
# print(c)
# print(random.sample(c,4)) #从列表a随机取b个数

str='i\'m singing in the rain'
print(str)
print(len(str)) #字符串长度
print(str.count('i')) #字符串中某个字符的个数
print(str.capitalize()) #首字母大写
print(str.isalpha())
print(str.isdigit())
print(str.islower())
print(str.istitle())
print(str.find('g')) #查找字符串中某个字符
print(str.replace('i','j')) #替换字符串中某个字符
l=str.split(' ') #用空格分割字符串，返回列表
print(l)
print(str.index('i')) #返回字符串中某字符的第一个索引
k=' '.join(l) #用空格拼接列表成为字符串
print(l,k)
print(type(l))
a1= [1,2,3]
l1=l.copy()
l1.append(a1) #列表末尾添加元素
l2=l.copy()
l2.extend(a1) #列表末尾拼接元素
print(l1)
print(l2)
l3=l.copy()
print(l3+a1) #同extend，但生成新的对象
l4=l3.insert(1,'not') #在索引位置处插入元素
print(l3)
print(type(k))
print(list(k)) #强制转换成列表
