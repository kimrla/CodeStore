import pymysql

# 打开数据库
# db=pymysql.connect(host='localhost',user='root',password='123456',database='tmp')
# #使用cursor方法创建游标对象
# cursor=db.cursor()
# #使用execute方法执行sql查询
# cursor.execute('select version()') #sql语句
# #使用fetch方法获取单条数据
# data=cursor.fetchone()
# print('database version:%s' % data)
# #关闭数据库
# db.close()

# # with语句连接数据库
# with pymysql.connect(host='localhost',user='root',password='123456',database='tmp') as db:
#     # 使用cursor方法创建游标对象
#     cursor = db.cursor()
#     # 使用execute方法执行sql查询
#     cursor.execute('select version()')  # sql语句
#     # 使用fetch方法获取单条数据
#     data = cursor.fetchone()
#     print('database version:%s' % data)

# # 创建数据库表
# with pymysql.connect(host='localhost',user='root',password='123456',database='tmp') as db:
#     # 使用cursor方法创建游标对象
#     cursor = db.cursor()
#     # 使用execute方法执行sql语句
#     cursor.execute('drop table if exists mytable')  # sql语句
#     # 预处理语句创建
#     sql='''create table mytable(
#         first_name char(20) not null,
#         last_name char(20) not null,
#         age int
#     )'''
#     cursor.execute(sql)

# # 插入数据
# with pymysql.connect(host='localhost', user='root', password='123456', database='tmp') as db:
#     # 使用cursor方法创建游标对象
#     cursor = db.cursor()
#     # 预处理语句创建
#     sql = '''insert into mytable(first_name,last_name,age)
#         values ('%s','%s',%s)''' % ('mac', 'mohan', 20)
#
#     try:
#         cursor.execute(sql)
#         db.commit()
#         print('插入成功')
#     except:
#         db.rollback()
#         print('插入失败，回滚')

# 查询数据
# with pymysql.connect(host='localhost', user='root', password='123456', database='tmp') as db:
#     # 使用cursor方法创建游标对象
#     cursor = db.cursor()
#     # 预处理语句创建
#     sql = '''select * from emp where sal > %s''' % 1000
#     try:
#         cursor.execute(sql)
#         results = cursor.fetchall()
#         print('查询结果共{}条'.format(cursor.rowcount))
#         for row in results:
#             print(row)
#     except:
#         print('查询发生错误')

# # 调用存储过程
# with pymysql.connect(host='localhost', user='root', password='123456', database='tmp') as db:
#     # 使用cursor方法创建游标对象
#     cursor = db.cursor()
#     cursor.callproc('get_maxsal')
#     result = cursor.fetchall()
#     print(cursor.rowcount)
#     print(result)

# 调用带参数的存储过程
with pymysql.connect(host='localhost', user='root', password='123456', database='tmp') as db:
    # 使用cursor方法创建游标对象
    cursor = db.cursor()
    #python只支持in参数，out和inout参数储存在服务器，用select查询，对应参数为@_储存过程名_0,1,2,3...
    cursor.callproc('get_sal',('BLAKE',0))
    cursor.execute('select @_get_sal_0,@_get_sal_1')
    print(cursor.fetchall())

