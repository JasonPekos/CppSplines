import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


dfSpline = pd.read_csv("output.csv") 
dfInput = pd.read_csv("input.csv")

plt.plot(dfSpline["t"],dfSpline["y"])
plt.plot(dfInput["t"], dfInput["y"], ".")


plt.show()


