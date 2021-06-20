from .base import Creator
from .. import param

from refractor import framework as rf

class DefaultConstants(Creator):

    def create(self):
        return rf.DefaultConstant()
