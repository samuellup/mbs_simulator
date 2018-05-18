import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-meta_in', action="store", dest='meta_in', required=True)
parser.add_argument('-out', action="store", dest='out', required=True)
args = parser.parse_args()

tag_list=list()

with open(args.meta_in, "r") as data:
	for line in data:
		if not line.startswith("#"):
			sp = line.split('\t')
			tag = sp[0] + "-" + sp[1]
			if tag not in tag_list:
				tag_list.append(tag)


sum_list=list()

for tag in tag_list:
	epsilon_rd=0
	epsilon_span=0
	n_samples=0
	with open(args.meta_in, "r") as data:
		for line in data:
			if not line.startswith("#"):
				sp2 = line.split('\t')
				tag2 = sp2[0] + "-" + sp2[1]
				if tag2 == tag:
					epsilon_rd = epsilon_rd + int(sp2[2])
					epsilon_span = epsilon_span + int(sp2[3])
					n_samples = n_samples + 1

	average_rd = float(epsilon_rd)/n_samples
	average_span = float(epsilon_span)/n_samples
	sub_list = [tag, average_rd, average_span]
	sum_list.append(sub_list)

with open(args.out, "w") as out:
	out.write("#RD" + "\t" + "MPS" + "\t" + "average_rd" + "\t" + "average_span" + "\n")

	for value in sum_list:
		sp = value[0].split("-")
		rd = sp[0]
		mps = sp[1]

		out.write(str(rd) + "\t" + str(mps) + "\t" + str(round(value[1], 1)) + "\t" + str(round(value[2], 1)) + "\n")