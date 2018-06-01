
from numpy.random import rand
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-in_va', action="store", dest='in_va', required=True)
parser.add_argument('-out', action="store", dest='out', required=True)
args = parser.parse_args()

fig, ax = plt.subplots()

#n = 25
#x, y = rand(2, n)
scale = 15

coord = list()
coords = list()
max_x=0
with open(args.in_va, "r") as variants:
	for line in variants:
		if not line.startswith("#"):
			coord = list()
			sp = line.split()
			position = int(sp[1])
			if position > max_x: max_x = position
			af = float(sp[6])/(float(sp[6]) + float(sp[5]))
			coord=[position, af]
			coords.append(coord)

for i in coords:
	ax.scatter(i[0], i[1], c='blue', s=scale, label='blue',
         alpha=0.3, edgecolors='blue')


#ax.legend()
ax.grid(True)

#plt.show()

plt.xlim(0, max_x)
plt.ylim(0, 1.1)

plt.savefig(args.out)
#plt.savefig('foo.pdf')

plt.close(fig)
