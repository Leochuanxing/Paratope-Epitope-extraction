'''
Analyze the middle aa
'''
import json
import os
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")
with open('middle_homo_Ab', 'r') as f:
    middle_homo_Ab = json.load(f)
with open('middle_homo_Ag', 'r') as f:
    middle_homo_Ag = json.load(f)
with open('middle_mouse_Ab', 'r') as f:
    middle_mouse_Ab = json.load(f)
with open('middle_mouse_Ag', 'r') as f:
    middle_mouse_Ag = json.load(f)

middle_mouse_Ag
middle_homo_Ag

middle_mouse_Ab
middle_homo_Ab

def Combine_middles(middle_homo, middle_mouse):
    middle = []
    for h in middle_homo:
        for i in middle_mouse:
            if h[0] == i[0]:
                h[1] += i[1]
        middle.append(h)
    return middle

middle_Ab = Combine_middles(middle_homo_Ab, middle_mouse_Ab)
middle_Ag = Combine_middles(middle_homo_Ag, middle_mouse_Ag)

middle_Ab.sort(key = lambda x:x[1])
middle_Ag.sort(key = lambda x:x[1])

middle_Ab
middle_Ag

len(middle_Ab)
len(middle_Ag)
aa20 = []
for i in middle_Ag:
    aa20.append(i[0])
len(aa20)  
middle_Ab.append(['CYS', 0])

n_middle_Ab = []
n_middle_Ag = []

for i in aa20:
    for j in middle_Ab:
        if i == j[0]:
            n_middle_Ab.append(j[1])
    for k in middle_Ag:
        if i == k[0]:
            n_middle_Ag.append(k[1])
aa20
n_middle_Ag            
middle_Ag
n_middle_Ab
middle_Ab
import matplotlib as plt
import numpy as np
x_pos =np.arange(20)

plt.figure(figsize=(20,10))
plt.bar(x_pos+0, n_middle_Ag, color='b', width = 0.3)
plt.bar(x_pos+0.3, n_middle_Ab, color='r', width = 0.3)
plt.xticks(x_pos, aa20)

#plt.set_size_inches(3, 1.5)
plt.show()

