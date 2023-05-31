import seaborn as sns
import matplotlib.pyplot as plt

penguins = sns.load_dataset("penguins")
sns.ecdfplot(data=penguins, x="flipper_length_mm")
sns.ecdfplot(data=penguins, y="flipper_length_mm")
sns.ecdfplot(data=penguins.filter(like="bill_", axis="columns"))

plt.figure()
sns.set_context("talk")
width = 6.5
height = 6.0
# width = 4
# height = 5

fig, ax = plt.subplots(figsize=(width, height))

ax = sns.ecdfplot(data=penguins, x="bill_length_mm", hue="species")
ax.legend(loc="center right", bbox_to_anchor=(1.4, 0.5))#,fancybox=True, ncol=1)
sns.despine()

plt.show()

