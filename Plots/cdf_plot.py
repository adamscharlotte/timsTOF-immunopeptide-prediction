import seaborn as sns
import matplotlib.pyplot as plt

penguins = sns.load_dataset("penguins")
# sns.ecdfplot(data=penguins, x="flipper_length_mm")
# sns.ecdfplot(data=penguins, y="flipper_length_mm")
# sns.ecdfplot(data=penguins.filter(like="bill_", axis="columns"))

sns.set_context("talk")
width = 6.5
height = 6.0
# width = 4
# height = 5
# penguins.species.drop_duplicates()
fig, ax = plt.subplots(figsize=(width, height))

ax = sns.ecdfplot(data=penguins, x="bill_length_mm", hue="species", palette=["#CDEAC0", "#f33b16", "#7A8DB3"])
ax.legend(loc="center right", bbox_to_anchor=(1.8, 0.5))#,fancybox=True, ncol=1)
sns.despine()
ax.grid(False, which="both")
plt.xlabel('Affinity score')

plt.show()

