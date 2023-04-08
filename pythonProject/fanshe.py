class dog:
    def __init__(self, name):
        self.name=name
    def eat(self):
        print('{} is eating'.format(self.name))
    def sleep(self):
        print(f'{self.name} is sleeping')

dog1=dog('tom')
# # dog1.eat()
# # dog1.sleep()
# cmd=input('请输入指令：')
# print(hasattr(dog1,cmd))
# if hasattr(dog1,cmd):
#     try:
#         getattr(dog1,cmd)()
#     except TypeError:
#         print(getattr(dog1,cmd))
# else:
#     print('%s can\'t do that'%dog1.name)
# print(dir(dog1))
# attr_key=input('请输入属性名：')
# attr_value=input('请输入属性值：')
# setattr(dog,attr_key,attr_value)
# print(dir(dog1))
# print(getattr(dog1,attr_key))

def say(content):
    print('tom说{}'.format(content))

method=input('请输入要绑定的方法名：')
setattr(dog1,method,say)
print(dir(dog1))
print(getattr(dog1,method)('一二三'))
delmet=input('请输入要删除的方法名：')
delattr(dog1,delmet)
print(dir(dog1))