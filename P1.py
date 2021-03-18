# %% Pilote 1

import scipy.io
from scipy.io import savemat
import glob
import os.path
import numpy as np


subs = ['sub-05DB', 'sub-07DB']

speeds = ['Normale', 'Rapide', 'Lente']

directions = ['Haut Bas', 'Bas Haut']

directions_order = ['Haut Bas', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Haut Bas', 'Bas Haut', 'Bas Haut', 'Bas Haut', 'Haut Bas', 'Haut Bas', 'Haut Bas', 'Bas Haut', 'Bas Haut']

speeds_order = ['Lente', 'Normale', 'Rapide', 'Lente', 'Normale', 'Lente', 'Normale', 'Rapide', 'Lente', 'Normale', 'Lente', 'Lente', 'Rapide', 'Normale', 'Rapide', 'Rapide', 'Normale', 'Lente', 'Rapide', 'Rapide', 'Lente', 'Rapide', 'Rapide', 'Rapide', 'Normale', 'Lente', 'Lente', 'Normale', 'Normale', 'Normale', 'Normale', 'Lente', 'Rapide', 'Normale', 'Lente', 'Normale', 'Lente', 'Lente', 'Normale', 'Lente', 'Normale', 'Normale', 'Rapide', 'Normale', 'Lente', 'Lente', 'Normale', 'Normale', 'Rapide', 'Normale', 'Lente', 'Rapide', 'Normale', 'Lente', 'Normale', 'Lente', 'Lente', 'Rapide', 'Normale', 'Normale', 'Lente', 'Rapide', 'Normale', 'Rapide', 'Lente', 'Normale', 'Lente', 'Rapide', 'Rapide', 'Rapide', 'Lente', 'Normale', 'Rapide', 'Rapide', 'Lente', 'Rapide', 'Lente', 'Lente', 'Lente', 'Rapide', 'Lente', 'Normale', 'Rapide', 'Rapide', 'Normale', 'Lente', 'Normale', 'Lente', 'Rapide', 'Lente', 'Normale', 'Normale', 'Normale', 'Normale', 'Normale', 'Rapide', 'Rapide', 'Rapide', 'Lente', 'Rapide', 'Lente', 'Normale', 'Rapide', 'Normale', 'Lente', 'Lente', 'Lente', 'Rapide', 'Rapide', 'Normale', 'Rapide', 'Rapide', 'Lente', 'Rapide', 'Lente', 'Rapide', 'Normale', 'Rapide', 'Rapide', 'Normale']

data=dict()


for sub in subs:
    
    idx = []
    datapath = '../Data/Données/Pilote 1 V2/'+ sub +'/'
    ##datapath = '../../Arthur/Pilote1/Données/'+ sub +'/'
    ##print(glob.glob(datapath + '*.mat'))
    for file in (glob.glob(datapath + '*.mat')):
        if os.path.basename(file)[-7]=='_':
            idx.append(os.path.basename(file)[-6:-4])
        else:

            idx.append(os.path.basename(file)[-7:-4])
        
    idx =  idx[:-14]
    idx=[int(x) for x in idx]
    idx.sort()
    for k in range(len(idx)):
        if idx[k]<100:
            idx[k]='%02d' % idx[k]
        else:
            idx[k]='%03d' % idx[k]
    ##print(sub)

    data[sub]=dict()
    for speed in speeds:
        data[sub][speed]=dict()
        ##print(speed)
        for direction in directions:
            data[sub][speed][direction]=dict()

    for k,l,m in zip(speeds_order, directions_order, idx):
        mat=scipy.io.loadmat(datapath + sub + '_' + m)
        data[sub][k][l][m]=mat['C3D'][0][0][2][0][0][2][:, 9:12]
        




# %%
import matplotlib.pyplot as plt

# plt.plot(data['Athina']['Lente']['Bas Haut']['47'][:, 0], color ='b', label ='x')
# plt.plot(data['Athina']['Lente']['Bas Haut']['47'][:, 1], color='r', label='y')

plt.plot(Vitesse['sub-07DB']['Lente']['Haut Bas']['67'], color='g', label='z')
plt.show()


#%% Butter filter

from scipy import signal
import math

order = 4
cutoff = 2

bb, ab = signal.butter(order, (2*cutoff)/100, 'low', analog=False, output='ba')


for sub in subs:
    for speed in speeds:
        for direction in directions:
            for idx in data[sub][speed][direction].keys():
                for k in range(3):

                    data[sub][speed][direction][idx][:,k] = signal.filtfilt(bb, ab, data[sub][speed][direction][idx][:,k])
                    #data[sub][speed][direction][idx][:,k] = signal.filtfilt(cb, ac, data[sub][speed][direction][idx][:,k])
                if  sub != 'Marie' and sub != 'Fatma' and sub != 'Athina' and sub != 'Elodie' and sub != 'Kevin':
                    data[sub][speed][direction][idx]=data[sub][speed][direction][idx][0:600]

# %%  del

# # sub-05DB
del data['sub-05DB']['Lente']['Haut Bas']['40']
del data['sub-05DB']['Lente']['Haut Bas']['57']
del data['sub-05DB']['Lente']['Haut Bas']['67']
del data['sub-05DB']['Lente']['Haut Bas']['86']
del data['sub-05DB']['Lente']['Haut Bas']['107']
del data['sub-05DB']['Lente']['Bas Haut']['71']
del data['sub-05DB']['Lente']['Bas Haut']['79']
del data['sub-05DB']['Lente']['Bas Haut']['115']
data['sub-05DB']['Lente']['Bas Haut']['27']=data['sub-05DB']['Lente']['Bas Haut']['27'][:500,:]

# sub-07DB
del data['sub-07DB']['Lente']['Haut Bas']['01']
del data['sub-07DB']['Lente']['Haut Bas']['09']
del data['sub-07DB']['Lente']['Haut Bas']['12']
del data['sub-07DB']['Lente']['Haut Bas']['18']
del data['sub-07DB']['Lente']['Haut Bas']['21']
del data['sub-07DB']['Lente']['Haut Bas']['32']
del data['sub-07DB']['Lente']['Haut Bas']['35']
del data['sub-07DB']['Lente']['Haut Bas']['37']
del data['sub-07DB']['Lente']['Haut Bas']['40']
del data['sub-07DB']['Lente']['Haut Bas']['45']

del data['sub-07DB']['Lente']['Haut Bas']['57']
del data['sub-07DB']['Lente']['Haut Bas']['61']
del data['sub-07DB']['Lente']['Haut Bas']['67']
del data['sub-07DB']['Lente']['Haut Bas']['86']
del data['sub-07DB']['Lente']['Bas Haut']['27']
del data['sub-07DB']['Lente']['Bas Haut']['38']
del data['sub-07DB']['Lente']['Bas Haut']['51']
del data['sub-07DB']['Lente']['Bas Haut']['54']
del data['sub-07DB']['Lente']['Bas Haut']['65']
del data['sub-07DB']['Lente']['Bas Haut']['71']
del data['sub-07DB']['Lente']['Bas Haut']['79']
del data['sub-07DB']['Lente']['Bas Haut']['101']
del data['sub-07DB']['Lente']['Bas Haut']['113']
# %% Calcul de la vitesse verticale pour chaque frame

import numpy as np

Vitesse = dict()
Vmax = dict()
VitesseQ = dict()
VQmax=dict()

for sub in data.keys():
    Vitesse[sub]=dict()
    Vmax[sub]=dict()
    VitesseQ[sub]=dict()
    VQmax[sub]=dict()
    for speed in data[sub].keys():
        Vitesse[sub][speed]=dict()
        Vmax[sub][speed]=dict()
        VitesseQ[sub][speed]=dict()
        VQmax[sub][speed]=dict()
        for direction in data[sub][speed].keys():
            Vitesse[sub][speed][direction]=dict()
            Vmax[sub][speed][direction]=dict()
            VitesseQ[sub][speed][direction]=dict()
            VQmax[sub][speed][direction]=dict()
            for idx in data[sub][speed][direction].keys():
                Vitesse[sub][speed][direction][idx]=np.zeros(len(data[sub][speed][direction][idx][:, 2]))
                VitesseQ[sub][speed][direction][idx]=np.zeros(len(data[sub][speed][direction][idx][:, 2]))

                for (k, v) in enumerate(data[sub][speed][direction][idx][:, 2]):
                    if k !=0:
                        Vitesse[sub][speed][direction][idx][k]=abs((data[sub][speed][direction][idx][k, 2] - data[sub][speed][direction][idx][k-1, 2])/0.01)
                        VitesseQ[sub][speed][direction][idx][k]=np.sqrt(pow((data[sub][speed][direction][idx][k, 2] - data[sub][speed][direction][idx][k-1, 2])/0.01, 2) + pow((data[sub][speed][direction][idx][k, 1] - data[sub][speed][direction][idx][k-1, 1])/0.01, 2) + pow((data[sub][speed][direction][idx][k, 0] - data[sub][speed][direction][idx][k-1, 0])/0.01, 2))

                Vmax[sub][speed][direction][idx]=np.nanmax(Vitesse[sub][speed][direction][idx])
                VQmax[sub][speed][direction][idx]=np.nanmax(VitesseQ[sub][speed][direction][idx])






# %% Visualisation de la vitesse verticale 

import matplotlib.pyplot as plt

plt.plot(Vitesse['sub-07DB']['Lente']['Haut Bas']['01'])
plt.show()

# for k in data[sub]['Normale']['Haut Bas'].keys():

#     plt.plot(Vitesse[sub]['Normale']['Haut Bas'][k], color='g', label='z')
#     plt.axhline(y=Vmax[sub]['Normale']['Haut Bas'][k])
#     plt.show()



# %% Création des 4 variables d'intérêt avec Vitesse verticale

dataT = dict()
VitesseT=dict()
VitesseQT=dict()
MS = dict()
ME=dict()
MDparam = 0.05
MD=dict()
TPV=dict()
SR=dict()
amp=dict()

for sub in data.keys():
    dataT[sub] = dict()
    VitesseT[sub]=dict()
    VitesseQT[sub]=dict()
    MS[sub] = dict()
    ME[sub]=dict()
    MD[sub]=dict()
    TPV[sub]=dict()
    SR[sub]=dict()
    amp[sub]=dict()

    for speed in data[sub].keys():
        dataT[sub][speed]=dict()
        VitesseT[sub][speed]=dict()
        VitesseQT[sub][speed]=dict()
        MS[sub][speed]=dict()
        ME[sub][speed]=dict()
        MD[sub][speed]=dict()
        TPV[sub][speed]=dict()
        SR[sub][speed]=dict()
        amp[sub][speed]=dict()

        for direction in data[sub][speed].keys():
            dataT[sub][speed][direction]=dict()
            VitesseT[sub][speed][direction]=dict()
            VitesseQT[sub][speed][direction]=dict()
            MS[sub][speed][direction]=dict()
            ME[sub][speed][direction]=dict()
            MD[sub][speed][direction]=dict()
            TPV[sub][speed][direction]=dict()
            SR[sub][speed][direction]=dict()
            amp[sub][speed][direction]=dict()

            for idx in data[sub][speed][direction].keys():
                
                
                s=0
                while (Vitesse[sub][speed][direction][idx][s] != Vmax[sub][speed][direction][idx]):
                    s=s+1
                t=s
                
                s=0
                    
                while (Vitesse[sub][speed][direction][idx][t-s]>= MDparam* Vmax[sub][speed][direction][idx]):
                    s = s+1
                MS[sub][speed][direction][idx]=t-s
                s=0
                print([sub, speed, direction, idx, Vmax[sub][speed][direction][idx]])
                
                while (Vitesse[sub][speed][direction][idx][t+s]>= MDparam* Vmax[sub][speed][direction][idx])  :
                    s = s+1
                ME[sub][speed][direction][idx]= t+s
                
                if (sub=='sub-02PC') and (speed=='Normale') and (direction=='Haut Bas') and (idx=='03'):
                    print(t)
                    print(MS[sub][speed][direction][idx])
                    print(ME[sub][speed][direction][idx])
                    print(ME[sub][speed][direction][idx]-MS[sub][speed][direction][idx])
                    print(t-MS[sub][speed][direction][idx])


                TPV[sub][speed][direction][idx]= (t - MS[sub][speed][direction][idx])/100


                SR[sub][speed][direction][idx]= (t  - MS[sub][speed][direction][idx])/(ME[sub][speed][direction][idx] - MS[sub][speed][direction][idx])
                VitesseT[sub][speed][direction][idx]= Vitesse[sub][speed][direction][idx][MS[sub][speed][direction][idx]:ME[sub][speed][direction][idx]]
                dataT[sub][speed][direction][idx]= data[sub][speed][direction][idx][MS[sub][speed][direction][idx]:ME[sub][speed][direction][idx]]
                MD[sub][speed][direction][idx]= (ME[sub][speed][direction][idx]- MS[sub][speed][direction][idx])/100
                amp[sub][speed][direction][idx]=abs(dataT[sub][speed][direction][idx][-1, 2]-dataT[sub][speed][direction][idx][0, 2])

#%% Percept test

idx = '17'

accT = np.zeros(len(VitesseT['sub-02PC']['Normale']['Haut Bas'][idx])-1)
for k in range(len(VitesseT['sub-02PC']['Normale']['Haut Bas'][idx])-1):
    accT[k]= VitesseT['sub-02PC']['Normale']['Haut Bas'][idx][k+1]-VitesseT['sub-02PC']['Normale']['Haut Bas'][idx][k]
    

# %% Save for percept 



percept = dict()
percept['data']=dataT['sub-02PC']['Normale']['Haut Bas']['03'][:,2]
savemat('coords.mat', percept)
#%% Check

import matplotlib.pyplot as plt

for speed in speeds :
    for direction in directions:
        for idx in data[sub][speed][direction].keys():
            print([sub, speed, direction, idx])
            plt.plot(Vitesse[sub][speed][direction][idx], color='g', label='z')
            plt.axhline(y=Vmax[sub][speed][direction][idx])
            plt.axvline(x=MS[sub][speed][direction][idx])
            plt.axvline(x=ME[sub][speed][direction][idx])
            plt.show()
# plt.axvline(x=MS[sub]['Lente']['Haut Bas']['01'])
# plt.axvline(x=ME[sub]['Lente']['Haut Bas']['01'])

# for sub in subs:
#     sub = 'sub-05DB'
#     for speed in data[sub].keys():
#         for direction in directions: 
#             for idx in data[sub][speed][direction].keys():
#                 print([sub, speed, direction, idx, Vmax[sub][speed][direction][idx]])
#                 plt.plot(Vitesse[sub][speed][direction][idx], color='g', label='z')
#                 plt.axhline(y=Vmax[sub][speed][direction][idx])
#                 plt.axvline(x=MS[sub][speed][direction][idx])
#                 plt.axvline(x=ME[sub][speed][direction][idx])
#                 plt.show()

#%% Excel

# sub-05DB
# del data['sub-05DB']['Lente']['Haut Bas']['40']
# del data['sub-05DB']['Lente']['Haut Bas']['57']
# del data['sub-05DB']['Lente']['Haut Bas']['67']
# del data['sub-05DB']['Lente']['Haut Bas']['86']
# del data['sub-05DB']['Lente']['Haut Bas']['107']
# del data['sub-05DB']['Lente']['Bas Haut']['71']
# del data['sub-05DB']['Lente']['Bas Haut']['79']
#del data['sub-05DB']['Lente']['Bas Haut']['115']


MD['sub-05DB']['Lente']['Haut Bas']['40']='NA'
SR['sub-05DB']['Lente']['Haut Bas']['40']='NA'
Vmax['sub-05DB']['Lente']['Haut Bas']['40']='NA'
amp['sub-05DB']['Lente']['Haut Bas']['40']='NA'

MD['sub-05DB']['Lente']['Haut Bas']['57']='NA'
SR['sub-05DB']['Lente']['Haut Bas']['57']='NA'
Vmax['sub-05DB']['Lente']['Haut Bas']['57']='NA'
amp['sub-05DB']['Lente']['Haut Bas']['57']='NA'

MD['sub-05DB']['Lente']['Haut Bas']['67']='NA'
Vmax['sub-05DB']['Lente']['Haut Bas']['67']='NA'
SR['sub-05DB']['Lente']['Haut Bas']['67']='NA'
amp['sub-05DB']['Lente']['Haut Bas']['67']='NA'

MD['sub-05DB']['Lente']['Haut Bas']['86']='NA'
Vmax['sub-05DB']['Lente']['Haut Bas']['86']='NA'
SR['sub-05DB']['Lente']['Haut Bas']['86']='NA'
amp['sub-05DB']['Lente']['Haut Bas']['86']='NA'

MD['sub-05DB']['Lente']['Haut Bas']['107']='NA'
Vmax['sub-05DB']['Lente']['Haut Bas']['107']='NA'
SR['sub-05DB']['Lente']['Haut Bas']['107']='NA'
amp['sub-05DB']['Lente']['Haut Bas']['107']='NA'

MD['sub-05DB']['Lente']['Bas Haut']['71']='NA'
Vmax['sub-05DB']['Lente']['Bas Haut']['71']='NA'
SR['sub-05DB']['Lente']['Bas Haut']['71']='NA'
amp['sub-05DB']['Lente']['Bas Haut']['71']='NA'


MD['sub-05DB']['Lente']['Bas Haut']['115']='NA'
Vmax['sub-05DB']['Lente']['Bas Haut']['115']='NA'
SR['sub-05DB']['Lente']['Bas Haut']['115']='NA'
amp['sub-05DB']['Lente']['Bas Haut']['115']='NA'

MD['sub-05DB']['Lente']['Bas Haut']['79']='NA'
Vmax['sub-05DB']['Lente']['Bas Haut']['79']='NA'
SR['sub-05DB']['Lente']['Bas Haut']['79']='NA'
amp['sub-05DB']['Lente']['Bas Haut']['79']='NA'

#sub-07DB

vars=[MD, Vmax, SR, amp]
for k in vars:

       k['sub-07DB']['Lente']['Haut Bas']['01']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['09']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['12']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['18']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['21']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['32']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['35']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['37']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['40']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['45']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['57']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['61']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['67']='NA'
       k['sub-07DB']['Lente']['Haut Bas']['86']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['27']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['38']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['51']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['54']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['65']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['71']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['79']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['101']='NA'
       k['sub-07DB']['Lente']['Bas Haut']['113']='NA'

import csv

subs = [ 'sub-05DB', 'sub-07DB']

speeds = [ 'Normale', 'Lente', 'Rapide']

directions = ['Haut Bas', 'Bas Haut']

table =np.zeros((len(data['sub-07DB']['Normale']['Haut Bas'])+1, 25))

for sub in subs:
    
    name = 'P1V2 - ' + sub  + '.csv'

    with open(name, 'w', newline ='') as file2:
        writer=csv.writer(file2)
        writer.writerow(['','MDDN','MDUN','MDDS','MDUS', 'MDDF', 'MDUF', 'VMDN', 'VMUN', 'VMDS', 'VMUS', 'VMDF', 'VMUF', 'TPVDN', 'TPVUN', 'TPVDS', 'TPVUS', 'TPVDF', 'TPVUF', 'AMPDN', 'AMPUN', 'AMPDS', 'AMPUS', 'AMPDF', 'AMPUF'])
        
        for k in range(len(MD['sub-07DB']['Normale']['Haut Bas'])):
                print(k)
                writer.writerow([ '', MD[sub]['Normale']['Haut Bas'][list(MD[sub]['Normale']['Haut Bas'].keys())[k]], MD[sub]['Normale']['Bas Haut'][list(MD[sub]['Normale']['Bas Haut'].keys())[k]], MD[sub]['Lente']['Haut Bas'][list(MD[sub]['Lente']['Haut Bas'].keys())[k]], MD[sub]['Lente']['Bas Haut'][list(MD[sub]['Lente']['Bas Haut'].keys())[k]], MD[sub]['Rapide']['Haut Bas'][list(MD[sub]['Rapide']['Haut Bas'].keys())[k]], MD[sub]['Rapide']['Bas Haut'][list(MD[sub]['Rapide']['Bas Haut'].keys())[k]], Vmax[sub]['Normale']['Haut Bas'][list(MD[sub]['Normale']['Haut Bas'].keys())[k]], Vmax[sub]['Normale']['Bas Haut'][list(MD[sub]['Normale']['Bas Haut'].keys())[k]], Vmax[sub]['Lente']['Haut Bas'][list(MD[sub]['Lente']['Haut Bas'].keys())[k]], Vmax[sub]['Lente']['Bas Haut'][list(MD[sub]['Lente']['Bas Haut'].keys())[k]], Vmax[sub]['Rapide']['Haut Bas'][list(MD[sub]['Rapide']['Haut Bas'].keys())[k]], Vmax[sub]['Rapide']['Bas Haut'][list(MD[sub]['Rapide']['Bas Haut'].keys())[k]], SR[sub]['Normale']['Haut Bas'][list(MD[sub]['Normale']['Haut Bas'].keys())[k]], SR[sub]['Normale']['Bas Haut'][list(MD[sub]['Normale']['Bas Haut'].keys())[k]], SR[sub]['Lente']['Haut Bas'][list(MD[sub]['Lente']['Haut Bas'].keys())[k]], SR[sub]['Lente']['Bas Haut'][list(MD[sub]['Lente']['Bas Haut'].keys())[k]], SR[sub]['Rapide']['Haut Bas'][list(MD[sub]['Rapide']['Haut Bas'].keys())[k]], SR[sub]['Rapide']['Bas Haut'][list(MD[sub]['Rapide']['Bas Haut'].keys())[k]], amp[sub]['Normale']['Haut Bas'][list(MD[sub]['Normale']['Haut Bas'].keys())[k]], amp[sub]['Normale']['Bas Haut'][list(MD[sub]['Normale']['Bas Haut'].keys())[k]], amp[sub]['Lente']['Haut Bas'][list(MD[sub]['Lente']['Haut Bas'].keys())[k]], amp[sub]['Lente']['Bas Haut'][list(MD[sub]['Lente']['Bas Haut'].keys())[k]], amp[sub]['Rapide']['Haut Bas'][list(MD[sub]['Rapide']['Haut Bas'].keys())[k]], amp[sub]['Rapide']['Bas Haut'][list(MD[sub]['Rapide']['Bas Haut'].keys())[k]]])




#%% Plot Vmax
import matplotlib.pyplot as plt
import numpy as np

error= [np.std([Vmax['sub-01CB']['Normale']['Haut Bas'][k] for k in Vmax['sub-01CB']['Normale']['Haut Bas']])/np.mean([Vmax['sub-01CB']['Normale']['Haut Bas'][k] for k in Vmax['sub-01CB']['Normale']['Haut Bas']]), np.std([Vmax['sub-02PC']['Normale']['Haut Bas'][k] for k in Vmax['sub-02PC']['Normale']['Haut Bas']])/np.mean([Vmax['sub-02PC']['Normale']['Haut Bas'][k] for k in Vmax['sub-02PC']['Normale']['Haut Bas']]), np.std([Vmax['sub-01CB']['Normale']['Bas Haut'][k] for k in Vmax['sub-01CB']['Normale']['Bas Haut']])/np.mean([Vmax['sub-01CB']['Normale']['Bas Haut'][k] for k in Vmax['sub-01CB']['Normale']['Bas Haut']]), np.std([Vmax['sub-02PC']['Normale']['Bas Haut'][k] for k in Vmax['sub-02PC']['Normale']['Bas Haut']])/np.mean([Vmax['sub-02PC']['Normale']['Bas Haut'][k] for k in Vmax['sub-02PC']['Normale']['Bas Haut']])]
labels=['D sub1', 'D sub 2', 'U sub 1', 'U sub 2']
# Build the plot
fig, ax = plt.subplots()
ax.bar([0, 1, 2, 3], [np.mean([Vmax['sub-01CB']['Normale']['Haut Bas'][k] for k in Vmax['sub-01CB']['Normale']['Haut Bas']]), np.mean([Vmax['sub-02PC']['Normale']['Haut Bas'][k] for k in Vmax['sub-02PC']['Normale']['Haut Bas']]), np.mean([Vmax['sub-01CB']['Normale']['Bas Haut'][k] for k in Vmax['sub-01CB']['Normale']['Bas Haut']]), np.mean([Vmax['sub-02PC']['Normale']['Bas Haut'][k] for k in Vmax['sub-02PC']['Normale']['Bas Haut']])], yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('Vmax')
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(labels)
ax.set_title('Pilote 1 : Vmax, Vitesse Normale')
ax.get_children()[3].set_color('orange')
ax.get_children()[4].set_color('orange')
plt.ylim([0,1200])
# ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
plt.savefig('pilot1Vmax.png')
plt.show()

#%% Plot MD
import numpy as np

error= [np.std([MD['sub-01CB']['Normale']['Haut Bas'][k] for k in MD['sub-01CB']['Normale']['Haut Bas']])/np.mean([MD['sub-01CB']['Normale']['Haut Bas'][k] for k in MD['sub-01CB']['Normale']['Haut Bas']]), np.std([MD['sub-02PC']['Normale']['Haut Bas'][k] for k in MD['sub-02PC']['Normale']['Haut Bas']])/np.mean([MD['sub-02PC']['Normale']['Haut Bas'][k] for k in MD['sub-02PC']['Normale']['Haut Bas']]), np.std([MD['sub-01CB']['Normale']['Bas Haut'][k] for k in MD['sub-01CB']['Normale']['Bas Haut']])/np.mean([MD['sub-01CB']['Normale']['Bas Haut'][k] for k in MD['sub-01CB']['Normale']['Bas Haut']]), np.std([MD['sub-02PC']['Normale']['Bas Haut'][k] for k in MD['sub-02PC']['Normale']['Bas Haut']])/np.mean([MD['sub-02PC']['Normale']['Bas Haut'][k] for k in MD['sub-02PC']['Normale']['Bas Haut']])]
labels=['D sub1', 'D sub 2', 'U sub 1', 'U sub 2']
# Build the plot
fig, ax = plt.subplots()
ax.bar([0, 1, 2, 3], [np.mean([MD['sub-01CB']['Normale']['Haut Bas'][k] for k in MD['sub-01CB']['Normale']['Haut Bas']]), np.mean([MD['sub-02PC']['Normale']['Haut Bas'][k] for k in MD['sub-02PC']['Normale']['Haut Bas']]), np.mean([MD['sub-01CB']['Normale']['Bas Haut'][k] for k in MD['sub-01CB']['Normale']['Bas Haut']]), np.mean([MD['sub-02PC']['Normale']['Bas Haut'][k] for k in MD['sub-02PC']['Normale']['Bas Haut']])], yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('MD')
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(labels)
ax.set_title('Pilote 1 : MD, Vitesse Normale')
ax.get_children()[3].set_color('orange')
ax.get_children()[4].set_color('orange')
plt.ylim([0,1.6])
# ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
plt.savefig('pilot1MD.png')
plt.show()

#%% Plot SR

import numpy as np

error= [np.std([SR['sub-01CB']['Normale']['Haut Bas'][k] for k in SR['sub-01CB']['Normale']['Haut Bas']])/np.mean([SR['sub-01CB']['Normale']['Haut Bas'][k] for k in SR['sub-01CB']['Normale']['Haut Bas']]), np.std([SR['sub-02PC']['Normale']['Haut Bas'][k] for k in SR['sub-02PC']['Normale']['Haut Bas']])/np.mean([SR['sub-02PC']['Normale']['Haut Bas'][k] for k in SR['sub-02PC']['Normale']['Haut Bas']]), np.std([SR['sub-01CB']['Normale']['Bas Haut'][k] for k in SR['sub-01CB']['Normale']['Bas Haut']])/np.mean([SR['sub-01CB']['Normale']['Bas Haut'][k] for k in SR['sub-01CB']['Normale']['Bas Haut']]), np.std([SR['sub-02PC']['Normale']['Bas Haut'][k] for k in SR['sub-02PC']['Normale']['Bas Haut']])/np.mean([SR['sub-02PC']['Normale']['Bas Haut'][k] for k in SR['sub-02PC']['Normale']['Bas Haut']])]
labels=['D sub1', 'D sub 2', 'U sub 1', 'U sub 2']
# Build the plot
fig, ax = plt.subplots()
ax.bar([0, 1, 2, 3], [np.mean([SR['sub-01CB']['Normale']['Haut Bas'][k] for k in SR['sub-01CB']['Normale']['Haut Bas']]), np.mean([SR['sub-02PC']['Normale']['Haut Bas'][k] for k in SR['sub-02PC']['Normale']['Haut Bas']]), np.mean([SR['sub-01CB']['Normale']['Bas Haut'][k] for k in SR['sub-01CB']['Normale']['Bas Haut']]), np.mean([SR['sub-02PC']['Normale']['Bas Haut'][k] for k in SR['sub-02PC']['Normale']['Bas Haut']])], yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('SR')
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(labels)
ax.set_title('Pilote 1 : TPV, Vitesse Normale')
ax.get_children()[3].set_color('orange')
ax.get_children()[4].set_color('orange')
plt.ylim([0,0.6])
# ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
plt.savefig('pilot1TPV.png')
plt.show()

#%% Plot Amp

import numpy as np

error= [np.std([amp['sub-01CB']['Normale']['Haut Bas'][k] for k in amp['sub-01CB']['Normale']['Haut Bas']])/np.mean([amp['sub-01CB']['Normale']['Haut Bas'][k] for k in amp['sub-01CB']['Normale']['Haut Bas']]), np.std([amp['sub-02PC']['Normale']['Haut Bas'][k] for k in amp['sub-02PC']['Normale']['Haut Bas']])/np.mean([amp['sub-02PC']['Normale']['Haut Bas'][k] for k in amp['sub-02PC']['Normale']['Haut Bas']]), np.std([amp['sub-01CB']['Normale']['Bas Haut'][k] for k in amp['sub-01CB']['Normale']['Bas Haut']])/np.mean([amp['sub-01CB']['Normale']['Bas Haut'][k] for k in amp['sub-01CB']['Normale']['Bas Haut']]), np.std([amp['sub-02PC']['Normale']['Bas Haut'][k] for k in amp['sub-02PC']['Normale']['Bas Haut']])/np.mean([amp['sub-02PC']['Normale']['Bas Haut'][k] for k in amp['sub-02PC']['Normale']['Bas Haut']])]
labels=['D sub1', 'D sub 2', 'U sub 1', 'U sub 2']
# Build the plot
fig, ax = plt.subplots()
ax.bar([0, 1, 2, 3], [np.mean([amp['sub-01CB']['Normale']['Haut Bas'][k] for k in amp['sub-01CB']['Normale']['Haut Bas']]), np.mean([amp['sub-02PC']['Normale']['Haut Bas'][k] for k in amp['sub-02PC']['Normale']['Haut Bas']]), np.mean([amp['sub-01CB']['Normale']['Bas Haut'][k] for k in amp['sub-01CB']['Normale']['Bas Haut']]), np.mean([amp['sub-02PC']['Normale']['Bas Haut'][k] for k in amp['sub-02PC']['Normale']['Bas Haut']])], yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('amp')
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(labels)
ax.set_title('Pilote 1 : amp, Vitesse Normale')
ax.get_children()[3].set_color('orange')
ax.get_children()[4].set_color('orange')
plt.ylim([0,700])
# ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
plt.savefig('pilot1amp.png')
plt.show()


# %%
