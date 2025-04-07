class Parent:
    def __init__(self, name):
        self.name = name
    def print_name(self, ChildA):
        print(self.name)
        print(self.sex)
class ChildA(Parent):
    def __init__(self, name, sex):
        super().__init__(name)
        # print_name(name)

    def do_something_with_b(self, b_instance):
        # 透過實例 b_instance 來呼叫 ChildB 的函式
        print(f"ChildA: 我可以呼叫 ChildB 的方法：{b_instance.some_method_in_b()}")

class ChildB(Parent):
    def __init__(self, name):
        super().__init__(name)

    def some_method_in_b(self):
        return "我是ChildB的方法"

# 建立實例並呼叫
a = ChildA("A的名字",'boy')
a.print_name(a)
# b = ChildB("B的名字")
# a.do_something_with_b(b)
