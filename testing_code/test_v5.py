class father():
    def __init__(self):
        self.eye = 2
        self.ear = 2
        self.nose = 1
        self.mouth = 1

class son(father):
    def __init__(self):   # 使用了 __init 的方法
        self.a = 100

class son2(father):
    def __init__(self):   # 使用了 __init 的方法
        super().__init__()
        print(self.eye)

oxxo = father()
print(oxxo.a)

