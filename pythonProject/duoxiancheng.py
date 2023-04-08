import threading #多线程模块
import time

namelist=['woniu','python','test']

def action(content):
    for item in content:
        print(threading.current_thread().getName()+item)#获取当前进程名+元素
        time.sleep(1)

thread=[]
for i in range(1,4):
    #target为要执行的函数名，args为函数的参数，参数必须以元组形式传入，若只有一个参数须用逗号
    t = threading.Thread(target=action,args=(namelist, ))
    t.daemon=True
    t.start() #启动线程
    thread.append(t)

# for t in thread:
#     t.join()
print('主线程结束')