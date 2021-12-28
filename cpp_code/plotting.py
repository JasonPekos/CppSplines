import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

dfSpline = pd.read_csv("cpp_code/output.csv")
dfInput  = pd.read_csv("cpp_code/input.csv")



plt.plot(dfSpline["t"], dfSpline["y"])
plt.plot(dfInput.iloc[:, 0], dfInput. iloc[:, 1], '.' )

plt.show()

