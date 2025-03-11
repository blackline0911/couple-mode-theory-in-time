class Parent:
    """父類別，所有子類別繼承自它"""
    pass

class ChildA(Parent):
    def __init__(self):
        self.attr_a = "我是子類 A"

    def update_child_b(self, other_instance):
        """使用 setattr() 修改另一個子類別的屬性"""
        if isinstance(other_instance, ChildB):  # 確保是 ChildB
            setattr(other_instance, "attr_b", "我是 A 修改的屬性")
        else:
            raise TypeError("只能修改 ChildB 的屬性")

class ChildB(Parent):
    def __init__(self):
        self.attr_b = "我是子類 B"

# 創建實例
a = ChildA()
b = ChildB()

print("修改前:", b.attr_b)  # ✅ "我是子類 B"

# 讓 ChildA 修改 ChildB 的屬性
a.update_child_b(b)

print("修改後:", b.attr_b)  # ✅ "我是 A 修改的屬性"
