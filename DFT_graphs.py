import matplotlib.pyplot as plt # importing relevant libraries for plotting
from matplotlib import colormaps
import numpy as np
import csv #to read csv file
import pandas as pd

#Reading files
#sig file
f_sig = open("original_sig.csv", "r") 
csv_reader = csv.reader(f_sig)
sig_data= pd.read_csv("original_sig.csv")
print('original signal data:', sig_data)
#dft file
f_dft = open("DFT_result.csv", "r")
csv_reader = csv.reader(f_dft)
dft_data= pd.read_csv("DFT_result.csv")
print('DFT data:', dft_data)


#Sorting data into dataframes
#sig file
df_sig = pd.DataFrame(sig_data)
#dft file
df_dft = pd.DataFrame(dft_data)

#statistics
print('dft statistics:', df_dft[['real_dft']].describe())

#Plotting
#3D plot
fig_3D = plt.figure()
ax = fig_3D.add_subplot(projection='3d')
#data
ax.scatter(df_sig.iloc[:, 0], df_sig.iloc[:, 2], df_sig.iloc[:, 1], s=3.14, color='pink', label='Original signal') # sig file
ax.plot(xs=df_dft.iloc[:, 0], ys=0,zs=df_dft.iloc[:, 1], color= 'deepskyblue', linewidth=1, label='DFT transformation') #dft file
#labels, title and legends
plt.title('Discrete Fourier Transform (DFT) Graph:', family='serif')
L = ax.legend(loc='upper right')
plt.setp(L.texts, family='serif')
plt.xlabel('N, #', family='serif')
plt.ylabel('imaginary, #', family='serif')
ax.set_zlabel('real, #', family='serif')

#show plt
plt.show()

#Closing files
f_sig.close() 
f_dft.close()