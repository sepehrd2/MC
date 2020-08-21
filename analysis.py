import function as fun
import numpy     as np
import math      as m


natom     = 64
lbox      = 4.0 
nstep     = 5000
nbins     = 100
dr        = 0.02
max_dis   = lbox * m.sqrt(3.0)
distances = np.zeros((natom * (natom - 1) * 0.5))
g         = np.zeros(nbins)
d         = np.zeros(3)
maxn      =  5
max_dis   = lbox * m.sqrt(3.0)
kvecs     = np.zeros(((maxn + 1)**3 , 3))
S_F       = np.zeros(len(kvecs))    
energy    = np.zeros(nstep)
T     = 3.0
INPUT1  = open("positions_"      + str(T) + ".txt" , "r")
INPUT2  = open("energy_"         + str(T) + ".txt" , "r")
OUTPUT5 = open("pair_cor_"       + str(T) + ".txt" , "w")
OUTPUT6 = open("stru_fac_"       + str(T) + ".txt" , "w")
data1   = INPUT1.read().split()
data2   = INPUT2.read().split()
r_last  = np.zeros((natom,3))

print '1- done with the reading'

for t in range(nstep - 1,nstep):
	for i in range(0,natom):
		r_last[i][0] = float(data1[t * natom * 3 + i * 3])
		r_last[i][1] = float(data1[t * natom * 3 + i * 3 + 1])
		r_last[i][2] = float(data1[t * natom * 3 + i * 3 + 2])

for t in range(0 ,nstep):
	energy[t] = float(data2[t * 2 + 1])
	print energy[t]

counter = 0

for i in range(0, natom):
	for j in range(i + 1, natom):
		d[0] = r_last[i][0] - r_last[j][0]
		d[1] = r_last[i][1] - r_last[j][1]
		d[2] = r_last[i][2] - r_last[j][2]

		d    = fun.my_disp_in_box(d, lbox)
		distances[counter] = fun.my_distance(d)

		if distances[counter] > max_dis:
			sys.exit("wrong distances :-(")

		counter = counter + 1
print '2- done with the distances'

g = fun.my_pair_correlation(distances, natom, nbins, dr, lbox)

print '3- done with the pair_cor'

for i in range(0, nbins):
	OUTPUT5.write('{0:.8f}  '.format(i * dr))
	OUTPUT5.write('{0:.8f}\n'.format(g[i]))

kvecs = fun.my_legal_kvecs(maxn, lbox)
S_F   = fun.my_calc_sk(kvecs, r_last)

print '4- done with the S_F'

for i in range(0, len(kvecs)):
	OUTPUT6.write('{}       '.format(i))
	OUTPUT6.write('{0:.8f}\n'.format(S_F[i]))
print "mean:   " + str(fun.my_mean(energy))
print "var:    " + str(fun.my_var(energy))
print "error:  " + str(fun.my_stderr(energy)) 
print "CL:     " + str(fun.my_actime(energy)) 

