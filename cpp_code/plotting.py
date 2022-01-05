import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set_theme() #Using Seaborn for aesthetics.
sns.set(rc={'figure.figsize': (10, 3)}) # Plot Dimensions.

sns.set_style("white") #Plot style.


dfSpline = pd.read_csv("output.csv") 
dfInput = pd.read_csv("input.csv")


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


