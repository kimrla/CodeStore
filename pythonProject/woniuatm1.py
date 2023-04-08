username = ['zhangsan','lisi','wangwu']
password = ['admin1','admin2','admin3']
def reg():
    print('请输入用户名：')
    newusername=input()
    if newusername in username:
        print('用户名已存在，注册失败')
        return False
    else:
        print('请输入密码：')
        newpassword = input()
        if len(newpassword)<6:
            print('密码长度小于6位，注册失败')
            return False
        else:
            username.append(newusername)
            password.append(newpassword)
            print('注册成功，请登录')
            return True
def login():
    print('请输入用户名：')
    loginusr=input()
    print('请输入密码：')
    loginpsw=input()
    if loginusr in username and loginpsw == password[username.index(loginusr)]:
        print('登录成功')
    else:
        print('登录失败')
        login()
if __name__=='__main__':
    if reg():
        login()
