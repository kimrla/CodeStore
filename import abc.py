import abc
class Client():
	def operation(self):
		pass
	
class Component(metaclass = abc.ABCMeta):#注意制作抽象类必须另外指定元类
	#operation的方法     
	@abc.abstractmethod
	def operation(self):
		pass

	@abc.abstractmethod
	def getChildren(self):
		pass

class Leaf(Component):
	def operation(self):
		pass
	def getChildren(self):
		pass
class Composite(Component):
	def operation(self):
		pass
	def getChildren(self):
		pass

Client.operation(Component)
Composite.Component()