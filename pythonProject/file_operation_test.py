# with open(r'zen.txt', mode='r', encoding='utf-8') as f:
#     # print(f.read())
#     # print(f.readline(5))
#     # print(f.readlines()[:3])
#     for line in f:
#         print(line)

with open(r'newtxt.txt',mode='w',encoding='utf-8') as nf:
    nf.write('这是第一行\n')
    nf.writelines(['这是第二行\n','这是第三行\n'])