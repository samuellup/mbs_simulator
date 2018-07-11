
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in_va', action="store", dest='in_va', required=True)
parser.add_argument('-out', action="store", dest='out', required=True)
parser.add_argument('-mps', action="store", dest='mps', required=True)
parser.add_argument('-rd', action="store", dest='rd', required=True)
args = parser.parse_args()

candidates_98=list()
candidates_95=list()

#Default values
n_candidates_95=0
n_candidates_98=0

cr_span_95=0
cr_span_98=0

with open(args.in_va, "r") as variants:
	for line in variants:
		if not line.startswith("#"):
			sp = line.split()
			position = int(sp[1])
			af = float(sp[6])/(float(sp[6]) + float(sp[5]))
			if af > 0.95:
				coord=[position, af]
				candidates_95.append(coord)

with open(args.in_va, "r") as variants:
	for line in variants:
		if not line.startswith("#"):
			sp = line.split()
			position = int(sp[1])
			af = float(sp[6])/(float(sp[6]) + float(sp[5]))
			if af > 0.98:
				coord=[position, af]
				candidates_98.append(coord)


pos_min_95=float("inf")
pos_max_95=0

pos_min_98=float("inf")
pos_max_98=0


# With curation
for var in candidates_95:
	try:
		d_95 = var[0] - ant_95
		if d_95 < 5845000:
			if var[0] < pos_min_95: pos_min_95 = var[0]
			if var[0] > pos_max_95: pos_max_95 = var[0]
			cr_span_95 = pos_max_95 - pos_min_95
			if cr_span_95 == 0: cr_span_95 = 1
		ant_95 = var[0]

	except:
		ant_95 = var[0]


for var in candidates_98:
	try:
		d_98 = var[0] - ant_98 
		if d_98 < 5845000:
			if var[0] < pos_min_98: pos_min_98 = var[0]
			if var[0] > pos_max_98: pos_max_98 = var[0]
			cr_span_98 = pos_max_98 - pos_min_98
			if cr_span_98 == 0: cr_span_98 = 1
		ant_98 = var[0]

	except:
		ant_98 = var[0]

n_candidates_95 = len(candidates_95)
n_candidates_98 = len(candidates_98)

#print n_candidates_98, cr_span_98
#print n_candidates_95, cr_span_95

with open(args.out, "a") as meta:
	meta.write(str(args.rd) + "\t" + str(args.mps) + "\t" + str(n_candidates_95) + "\t" + str(cr_span_95) + "\t" + str(n_candidates_98) + "\t" + str(cr_span_98) + "\n")


########################################## Deprecated ######################################################################################################
'''
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
			if af > 0.95:
				coord=[position, af]
				candidates.append(coord)

pos_min=float("inf")
pos_max=0



# Without curation
'''
'''
for var in candidates:
	if var[0] < pos_min: pos_min = var[0]
	if var[0] > pos_max: pos_max = var[0]
	cr_span = pos_max - pos_min
	if cr_span == 0: cr_span = 1
'''
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

'''

########################################## Deprecated ######################################################################################################
