import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-nbr_mutations', action="store", dest = 'nbr_in')
args = parser.parse_args()



n_in = int(args.nbr_in)

# 1)  48-52 percentage values
low = int(0.48*n_in)
high = int(0.52*n_in)

# 2) Random number choice
from random import randint
rand = int(randint(low, high))

print rand




'''



nbr_mutations_1=`python2 an_scripts/rand.py -nbr_mutations $nbr_mutations 2>> $my_log_file`
nbr_mutations_2= $(( nbr_mutations-nbr_mutations_1 ))





'''