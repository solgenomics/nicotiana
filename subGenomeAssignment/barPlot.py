import numpy as np
import matplotlib.pyplot as plt

N = 10
Nsyl = np.array([797, 2523, 1785, 2379, 3749, 1669, 2289, 869, 3119, 4176])
Ntom = np.array([919, 3038, 2235, 2775, 4314, 1922, 2667, 1045, 3687, 4893])
#insig = np.array([12885, 12787, 13212, 12384, 12403, 11657, 11655, 11577, 12662, 12319]) - Nsyl - Ntom

xpos = np.arange(N)
width = 0.25
p1 = plt.bar(xpos, Nsyl, width, color='#0000CD')
p2 = plt.bar(xpos+width, Ntom, width, color='#FF0000')
#p3 = plt.bar(xpos, insig, width, bottom = (Nsyl+Ntom), color='#A9A9A9')
plt.ylabel('# of gene pairs')
plt.title('HEB for various tissue types')
plt.xticks(xpos, ('ShootApex','Shoot', 'Root', 'Juvenile_Leaf', 'Lower_Leaf', 'Mid_Leaf', 'Upper_Leaf', 'ImmatureFlowers', 'Petal', 'Sepal'))
plt.legend((p1[0], p2[0]), ('Nsyl', 'Ntom'), fontsize='x-large')

# add percentage for Nsyl
for bar in p1:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2.0, height/2, '%d' % int(height), ha='center', va='bottom', color='#FFFFFF')

# add percentage for Ntom
for i in range(len(p2)):
    bar1 = p1[i]
    bar2 = p2[i]
    height = p2[i].get_height()
    plt.text(bar1.get_x()+1.5*bar1.get_width(), height/2, '%d '% int(height), ha='center', va='bottom', color="#FFFFFF")
# add percentage for Non-sig


plt.show()
