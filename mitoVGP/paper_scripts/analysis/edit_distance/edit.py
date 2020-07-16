import random
from sys import argv
dna = set(["A", "C", "G", "T"])

def deletion(seq):
	chars = list(seq)
	idx = random.choice(range(len(chars)))
	del chars[idx]
	return "".join(chars)

def insertion(seq):
	chars = list(seq)
	idx = random.choice(range(len(chars)))
	new_base = random.choice(list(dna))
	chars.insert(idx, new_base)
	return "".join(chars)

def substitute(seq):
	chars = list(seq)
	idx = random.choice(range(len(chars)))
	new_base = random.choice(list(dna.difference(chars[idx])))
	chars[idx] = new_base
	return "".join(chars)

class Sequence(str):

	def mutate(self, d):
		mutant = self
		m = 0
		while m < d:
			if m == d:
				break
			mutant_type = str(random.choices(["d", "s", "i"], weights=[0.29, 0.29, 0.42],k=1)[0])
			if mutant_type == "i":
				mutant = insertion(mutant)
			elif mutant_type == "d":
				mutant = deletion(mutant)
			elif mutant_type == "s":
				mutant = substitute(mutant)
			m += 1
		return str(mutant)

#ref
arg1 = argv[1]
arg2 = argv[2]
arg3 = argv[3]
arg4 = argv[4]

with open(arg1, 'r') as f:
	lines = f.readlines()

	hd = lines[0].strip()
	rd = "".join(line.strip() for line in lines[1:])

s = Sequence(str(rd))
d = 1

for n in range(1, int(arg2)+1):

	fn = "%s/%s_%s.fasta" % (arg3, hd.replace('>', ''), d*int(arg4))
	f = open(fn, "w")
	edited = s.mutate(d*int(arg4))
	hd_new = "%s_%s" % (hd, d*int(arg4))
	f.write(hd_new+"\n")
	f.write(str(edited)+"\n")
	d += 1

f.close()
