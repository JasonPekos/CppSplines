import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set_theme()
sns.set(rc={'figure.figsize': (10, 3)})

sns.set_style("white")


dfSpline = pd.read_csv("cpp_code/output.csv")
dfInput = pd.read_csv("cpp_code/input.csv")


sns.lineplot(
    data=dfSpline,
    x=dfSpline["t"],
    y=dfSpline["y"],
    palette="crest")

sns.scatterplot(
    data=dfInput,
    x=dfInput["t"],
    y=dfInput["y"],
    palette="deep"
)


plt.show()
