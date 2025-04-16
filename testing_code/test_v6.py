class Parent:
    pass  # 這個父類不用做任何事

class ChildA(Parent):
    def __init__(self):
        self.attr_a = "我是子類A"

    def share_with(self, other_instance):
        """將自己的屬性傳給另一個子類別"""
        for key, value in self.__dict__.items():
            setattr(other_instance, key, value)

class ChildB(Parent):
    def __init__(self):
        self.attr_b = "我是子類B"

# 創建子類實例
a = ChildA()
b = ChildB()

# 讓 ChildA 傳遞屬性給 ChildB
a.share_with(b)

print(b.attr_a)  # ✅ "我是子類A" (成功獲取 ChildA 的屬性)
print(b.attr_b)  # ✅ "我是子類B"
