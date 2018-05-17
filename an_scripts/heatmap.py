import numpy as np
import numpy.random
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('-in_data', action="store", dest='in_data', required=True)
parser.add_argument('-out', action="store", dest='out', required=True)
parser.add_argument('-mps', action="store", dest='mps', required=True)
parser.add_argument('-rd', action="store", dest='rd', required=True)
args = parser.parse_args()



# Create heatmap
heatmap, xedges, yedges = np.histogram2d(x, y, bins=(mps,rd))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
 
# Plot heatmap
plt.clf()
plt.title('Pythonspot.com heatmap example')
plt.ylabel('y')
plt.xlabel('x')
plt.imshow(heatmap, extent=extent)
plt.show()