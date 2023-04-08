# def foo():
#     print('in foo')
#
# def boo(func):
#     print(func)
#     func()
#
# boo(foo)

# def foo():
#     sum=0
#     for i in range(10000000):
#         sum += i
#     print(sum)
# import  time
# def gf(func):
#     start_time=time.time()
#     func()
#     end_time=time.time()
#     print(f'运行时间为{end_time-start_time}秒')
#     print('运行时间为{}秒'.format(end_time-start_time))
#     print('运行时间为%.16f秒'%(end_time-start_time))
#
# gf(foo)

# def foo():
#     sum=0
#     for i in range(10000000):
#         sum += i
#     print(sum)
# import  time
# def gf(func):
#     start_time=time.time()
#     func()
#     end_time=time.time()
#     print(f'运行时间为{end_time-start_time}秒')
#     print('运行时间为{}秒'.format(end_time-start_time))
#     print('运行时间为%.16f秒'%(end_time-start_time))
#     return func
#
# foo=gf(foo)
# foo()

# def foo():
#     print('in foo')
#     def boo():
#         print('in boo')
#
# foo()
# def boo():
#     print('in boo')
#
# def foo(func):
#
#     def gf():
#         print('in foo')
#         func()
#     return gf
# boo=foo(boo)
# boo()
def foo(func):

    def gf():
        print('in foo')
        func()
    return gf
@foo
def boo():
    print('in boo')

boo()