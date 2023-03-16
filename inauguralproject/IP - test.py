#importing packages
from types import SimpleNamespace

import numpy as np
from scipy import optimize

import pandas as pd
import matplotlib.pyplot as plt

#Defining the class
class household : 
    def __init__(self):
        
        #creating namespaces
        par = self.par = SimpleNamespace()
        sol = self.sol = SimpleNamespace()
        