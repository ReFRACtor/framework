from .base import Creator
from .. import param

import refractor.framework as rf

class DefaultConstants(Creator):

    def create(self):
        return rf.DefaultConstant()
