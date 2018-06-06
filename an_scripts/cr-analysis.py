import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in_va', action="store", dest='in_va', required=True)
parser.add_argument('-out', action="store", dest='out', required=True)
parser.add_argument('-mps', action="store", dest='mps', required=True)
parser.add_argument('-rd', action="store", dest='rd', required=True)
args = parser.parse_args()

candidates=list()

#Default values
n_candidates=0
cr_span=0

with open(args.in_va, "r") as variants:
	for line in variants:
		if not line.startswith("#"):
			sp = line.split()
			position = int(sp[1])
			af = float(sp[6])/(float(sp[6]) + float(sp[5]))
			if af > 0.98:
				coord=[position, af]
				candidates.append(coord)

pos_min=float("inf")
pos_max=0



# Without curation
'''
for var in candidates:
	if var[0] < pos_min: pos_min = var[0]
	if var[0] > pos_max: pos_max = var[0]
	cr_span = pos_max - pos_min
	if cr_span == 0: cr_span = 1
'''

# With curation
for var in candidates:
	try:
		d = var[0] - ant 
		if d < 4000000:
			if var[0] < pos_min: pos_min = var[0]
			if var[0] > pos_max: pos_max = var[0]
			cr_span = pos_max - pos_min
			if cr_span == 0: cr_span = 1
		ant = var[0]

	except:
		ant = var[0]


n_candidates = len(candidates)

with open(args.out, "a") as meta:
	meta.write(str(args.rd) + "\t" + str(args.mps) + "\t" + str(n_candidates) + "\t" + str(cr_span) + "\n")


#print args.rd, args.mps, n_candidates, cr_span
