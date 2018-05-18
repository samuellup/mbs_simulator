import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in_data', action="store", dest='in_data', required=True)
parser.add_argument('-out', action="store", dest='out', required=True)
args = parser.parse_args()

# Read in the relevant data from our input file
dt = np.dtype([('rd', np.int), ('mps', np.int), ('av', np.float)])
data = np.genfromtxt('yo.txt', dtype=dt, usecols=(1,2,3),
                     delimiter=('\t'))
print data

'''
heatmap = np.empty((3, 3))

for rd, mps, av in data:
    # NumPy arrays are zero-indexed; days and months are not!
    heatmap[rd, mps] = av



# Plot the heatmap, customize and label the ticks
fig = plt.figure()
ax = fig.add_subplot(111)
im = ax.imshow(heatmap, interpolation='nearest')
#ax.set_yticks(range(12))
#ax.set_yticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
#                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
#days = np.array(range(0, 31, 2))
#ax.set_xticks(days)
#ax.set_xticklabels(['{:d}'.format(day+1) for day in days])
#ax.set_xlabel('Day of month')
#ax.set_title('Maximum daily temperatures in Boston, 2012')

# Add a colour bar along the bottom and label it
cbar = fig.colorbar(ax=ax, mappable=im, orientation='horizontal')
cbar.set_label('Aver')

plt.savefig(args.out)
'''