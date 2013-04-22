'''
Created on 18/04/2013

@author: JuanCarlos
'''

class CircularStack(object):
    '''
    classdocs
    '''
    data = []
    size = 0
    index = 0
    tail = 0


    def __init__(self,size = 0):
        '''
        Constructor
        '''
        if size > 0:
            self.data = [0] * size
        self.index = 0
        self.tail = 0
    
    def push(self,data):
        next_tail = self._get_next_tail()
        if next_tail == self.index:
            raise OverflowError
        else:
            self.data[self.tail] = data
            self.tail += 1

    def pop(self):
        if self.index == self.tail:
            raise IndexError
        else:
            result = self.data[self.index]
            self.index = (self.index + 1) % self.size
            return result
        
    def _get_next_tail(self):
        return (self.tail + 1) % self.size
        