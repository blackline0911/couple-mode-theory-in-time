class Animal():
    def __init__(self,name):
        self.class_name = name
        self.eyes_num = 2
        self.legs_num = 4
    def breath(self, animal):
        print(f"{animal}呼吸")
class Dog(Animal):
    def __init__(self):
        super().__init__(name="動物")
        self.class_name = "小狗"
    def breath(self, dog="小狗"):
        print(super().eyes_num)
        
dog = Dog()
print(f"{dog.class_name}有{dog.eyes_num}個眼睛、{dog.legs_num}隻腳")
dog.breath()