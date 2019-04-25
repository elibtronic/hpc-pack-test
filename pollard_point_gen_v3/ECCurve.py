
from collections import namedtuple
import random
import operator
import re
import sys
import itertools

from random import randint
from random import getrandbits
from random import shuffle
from random import seed
from math import sqrt
#import numpy as np
#from scipy.stats import chisquare



#Config Variables
## Evolution Parameters
POPSIZE = 500
NUMGEN = 100
CXPB = 0.9
MUTPB = 0.1
MINTREE = 1
MAXTREE = 2
TOURSIZE = 5
MAXTREESIZE = 17
MUT_TREE_MIN = 0
MUT_TREE_MAX = 2


# fn,a,b,n,P.x,P.y,Q.x,Q.y
#TC = [16077749,87,30,16072891,8947224,7340826,15387666,13125775]
#Yes this is ugly
try:
	TC = open(sys.argv[2],"r").readline().split(",")
	TC = [int(i) for i in TC]
except:
	print("No Curve File Specified")
	sys.exit()


## Log / Seed Parameters

ORIGINAL_FILE = "./data_prep/"+str(TC[0])+"_original_seed_"

POINTS_FILE_X = "./data_prep/"+str(TC[0])+"_x_points_seed_"
POINTS_FILE_Y = "./data_prep/"+str(TC[0])+"_y_points_seed_"

LONGRUN_FILE = "./data_prep/"+str(TC[0])+"_long_runs.txt"

MUL_A_FILE = "./data_prep/"+str(TC[0])+"_a_mul_seed_"
MUL_B_FILE = "./data_prep/"+str(TC[0])+"_b_mul_seed_"

C_PRIME_FILE = "./data_prep/"+str(TC[0])+"_cprime.txt"
D_PRIME_FILE = "./data_prep/"+str(TC[0])+"_dprime.txt"

GP_LOG_BASE = "./GP/"
GA_LOG_BASE = "./GA/"

GP_EVO_CACHE_BASE = "gp_evo_cache.txt"


CURVE_DEBUG = True
RHO_DEBUG = False
FIT_DEBUG = False
UBER_SEED = 23232

GP_SHOW_EVOL = False
GA_SHOW_EVOL = False


#Also known as L
PARTITIONS = 32

#GP Specific
NUM_TEST_POINTS = 512
PEN_THRESH = 17


#FORCE_PARTITION_FUNCTION = "(P.x)"
FORCE_PARTITION_FUNCTION =  'operator.mul(operator.add(P.y, 3), operator.mul(1, P.y))'

# Maximum multiple to consider for Random Point Generation
MAXPOINTRANGE = 1000000


#Rho mul and inits, cprime and dprime as set to interval

INIT_C = int(TC[3]*0.25)
INIT_D = int(TC[3]*0.75)


try:
	TEST_POINTS_X = [int(line.strip()) for line in open(POINTS_FILE_X+str(UBER_SEED)+".txt",'r')]
	TEST_POINTS_Y = [int(line.strip()) for line in open(POINTS_FILE_Y+str(UBER_SEED)+".txt",'r')]
except:
	pass



 
from collections import namedtuple
import random
import operator
import re
import sys
import itertools

from random import randint
from random import getrandbits
from random import shuffle
from random import seed
from math import sqrt
#import numpy as np
#from scipy.stats import chisquare


#from trigger_complete import *

#Point on the curve
Point = namedtuple("Point","x y")


## For calculating the actual value of Rho
#  Inverting group arithemtic
#Via https://rosettacode.org/wiki/Modular_inverse#Python
def extended_gcd(aa, bb):
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(lastremainder, remainder)
        x, lastx = lastx - quotient*x, x
        y, lasty = lasty - quotient*y, y
    return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)
def modinv(a, m):
	g, x, y = extended_gcd(a, m)
	if g != 1:
		raise ValueError
	return x % m


class ECCurve:
	""" Basic ECCurve with some elementary operations, including Rho
		F(_prime) y^2 = x^3 + a*x + b
	"""

	# Basic Algebra stuff
	Origin = 'Origin'
	prime = 0			# Field Size
	a = 0				# a param of Weisterass eq
	b = 0				# b param of Weisterass eq
	points = []			# List of the points, not calculated until .calc_all_points()
	cardinality = 0		# Or prime order of field, no built in methods to calc it yet written
						# ideally same value of n passed to Rho
						
	#Debugging only
	debug = True
	
	
	## Rho Information
	#
	#
	#When re-doing expected runs
	init_a = [] 		# P.X values for the initial seed points, must be self.L long
	init_b = []			# P.y values for the initial seed points, must be self.L long
	
	init_c_prime = -1	# Initial guess for c' in c'P + d'Q on first run
	init_d_prime = -1	# Initial guess for d' in c'P + d'Q on first run

	birthday_cutoff = 0		# Stop rho after the value the birthday paradox is reached

	random_points = []	#used by fitness
	
	P = False	# P point
	Q = False	# Q point

	R = []	#initial points vector constructed with a_i, b_i

	#Needed for Rho calculations
	partition_fn = ""	# Partition Function, calculated and then Mod self.L is always applied
	rho_iters = 0		# How many times it run
	L = 0				# How many subsections Rho will use
	k = 0				# Q=kP, won't be calculated unless specified

	rho_runs = 0

	point_buffer = []

	#GA Implementation
	partition_pattern = []	#Represent partition function as a list, where each index is the next partition to pick, 0 indexed
	pattern_next = 0	#Which position in the pattern to pick next 


	#For evolutionary stuff
	fitness_fn = ""
	fitness = 0.0
	fit_penalty = 0.0

	def __init__(self, prime,a,b):
		self.prime = prime
		self.a = a
		self.b = b

	def generate_init_muls(self,n):
		for j in range(0,self.L):
				self.init_a.append(randint(0,n-1))
		for j in range(0,self.L):
				self.init_b.append(randint(0,n-1))

	def generate_random_points(self,numPoints=0):
		"""
		Pre-Computes the points needed for fitness
		"""
		
		if numPoints == 0:
			up_to = self.L
		else:
			upt_to = numPoints
		
		for j in range(0,numPoints):
			self.random_points.append(self.rand_point())
		
		return self.random_points
		
	
	def load_l_points(self,numPoints,sortPoints=True):
		"""
		Loads first L points in points file
		"""
		pf = open(POINTS_FILE,"r")
		for i in range(0,numPoints):
			p = pf.readline()
			x_coord = int(p.split(",")[0])
			y_coord = int(p.split(",")[1].strip("\n"))
			self.random_points.append(Point(x_coord,y_coord))
		
		if sortPoints:
			self.random_points.sort()
		
		return self.random_points
		
	def load_mul(self):
		"""
		Depreciated: Loads pre-done a,b multipliers
		"""
		af = open(MUL_A_FILE,"r")
		bf = open(MUL_B_FILE,"r")
		
		for a in af.readlines():
			self.init_a.append(int(a))
			
		for b in bf.readlines():
			self.init_b.append(int(b))
		
		return True
		
	def load_random_points(self,numPoints):
		"""
		Decpreciated: Loads random points from file outlined in para file
		"""
		big_p_list = []
		pf = open(POINTS_FILE,"r")
		for p in pf.readlines():
			x_coord = int(p.split(",")[0])
			y_coord = int(p.split(",")[1].strip("\n"))
			big_p_list.append(Point(x_coord,y_coord))
		
		for i in range(0,numPoints):
			rp = random.choice(big_p_list)
			self.random_points.append(rp)
			big_p_list.remove(rp)
		
		#print(self.random_points)
		return self.random_points
	

		

	def rho(self,n,value=False,iter_limit=0):
		"""
		Serial Rho Implementation, will return k if last para is passed as True
		Set value to true to get the final value of the multiplier
		Set iter_limit to some number for how many times each rho run will iterate
		Uses GP representation
		"""
		#R = []
		a = []
		b = []
		
		self.rho_runs += 1
		
		a = self.init_a
		b = self.init_b
		c_prime = self.init_c_prime
		d_prime = self.init_d_prime

		#print(c_prime)
		#print(d_prime)

		if len(self.R)== 0:
			for j in range(0,self.L):
				self.R.append(self.add(self.mul(self.P,self.init_a[j]),self.mul(self.Q,self.init_b[j])))

		X_prime = self.add(self.mul(self.P,c_prime),self.mul(self.Q,d_prime))
		X_dbl_prime = X_prime
		c_dbl_prime = c_prime
		d_dbl_prime = d_prime
		
		
		if RHO_DEBUG == True:
			print("a"+str(self.init_a))
			print("b"+str(self.init_b))
			print("itr,c',d',X',c'',d'',X''")
		
		#Do I want this here?
		self.rho_iters = 0
		
		while True:
			
			self.rho_iters += 1
			
			#if self.rho_iters % 500 == 0:
			#	print(".", end="")
			
			#If we want to stop at Birthday Limit
			if self.birthday_cutoff != 0 and self.rho_iters > self.birthday_cutoff:
				self.rho_iters = False
				return False
			
			#For when we only want to 
			if iter_limit != 0:
				if self.rho_iters == (self.rho_runs * iter_limit):
					return False
			
			j = self.partition(X_prime)
			X_prime = self.add(X_prime,self.R[j])
			c_prime = (c_prime + a[j]) % n
			d_prime = (d_prime + b[j]) % n
			
			for i in range(1,2+1):
				j = self.partition(X_dbl_prime)
				X_dbl_prime = self.add(X_dbl_prime,self.R[j])
				c_dbl_prime = (c_dbl_prime + a[j]) % n
				d_dbl_prime = (d_dbl_prime + b[j]) % n 
			
			if RHO_DEBUG == True:
				print(str(self.rho_iters)+"," + str(c_prime) + ","+ str(d_prime) + "," + str(X_prime) + "," + str(c_dbl_prime) + ","+ str(d_dbl_prime) + "," + str(X_dbl_prime))
			
			if self.points_equal(X_prime,X_dbl_prime):
				break
		
		if d_prime == d_dbl_prime:
			return False
		else:
			if value == True:
				l = (c_prime - c_dbl_prime) * modinv((d_dbl_prime - d_prime),n) % n
				self.k = l
				return l
			else:
				return True

	def rho_ga(self,n,value=False,iter_limit=0):
		"""
		Serial Rho Implementation, will return k if last para is passed as True
		Set value to true to get the final value of the multiplier
		Set iter_limit to some number for how many times each rho run will iterate
		Uses GA representation of partition
		"""
		R = []
		a = []
		b = []
		
		self.rho_runs += 1
		
		a = self.init_a
		b = self.init_b
		c_prime = self.init_c_prime
		d_prime = self.init_d_prime

		for j in range(0,self.L):
			R.append(self.add(self.mul(self.P,self.init_a[j]),self.mul(self.Q,self.init_b[j])))

		X_prime = self.add(self.mul(self.P,c_prime),self.mul(self.Q,d_prime))
		X_dbl_prime = X_prime
		c_dbl_prime = c_prime
		d_dbl_prime = d_prime
		
		
		if RHO_DEBUG == True:
			print("itr,c',d',X',c'',d'',X''")
		
		#Do I want this here?
		self.rho_iters = 0
		
		while True:
			
			self.rho_iters += 1
			
			#if self.rho_iters % 500 == 0:
			#	print(".", end="")
			
			#If we want to stop at Birthday Limit
			if self.birthday_cutoff != 0 and self.rho_iters > self.birthday_cutoff:
				self.rho_iters = False
				return False
			
			#For when we only want to 
			if iter_limit != 0:
				if self.rho_iters == (self.rho_runs * iter_limit):
					return False
			
			j = self.partition_next(X_prime)
			X_prime = self.add(X_prime,R[j])
			c_prime = (c_prime + a[j]) % n
			d_prime = (d_prime + b[j]) % n
			
			for i in range(1,2+1):
				j = self.partition_next(X_dbl_prime)
				X_dbl_prime = self.add(X_dbl_prime,R[j])
				c_dbl_prime = (c_dbl_prime + a[j]) % n
				d_dbl_prime = (d_dbl_prime + b[j]) % n 
			
			if RHO_DEBUG == True:
				print(str(self.rho_iters)+"," + str(c_prime) + ","+ str(d_prime) + "," + str(X_prime) + "," + str(c_dbl_prime) + ","+ str(d_dbl_prime) + "," + str(X_dbl_prime))
			
			if self.points_equal(X_prime,X_dbl_prime):
				break
		
		if d_prime == d_dbl_prime:
			return False
		else:
			if value == True:
				l = (c_prime - c_dbl_prime) * modinv((d_dbl_prime - d_prime),n) % n
				self.k = l
				return l
			else:
				return True

	def build_R(self):
		for j in range(0,self.L):
			self.R.append(self.add(self.mul(self.P,self.init_a[j]),self.mul(self.Q,self.init_b[j])))


### OLD GP Fitness Functions
### Very vestigial

	#def calc_fitness_oo(self):
		#"""
		#Only want 1 per each partition, fit_penalty applied to each non-unique value
		#"""
		#scores = []
		#for P in self.random_points:
			#scores.append(self.partition(P))
		#missing_penalty = float((self.L - len(set(scores))) * self.fit_penalty)
		#self.fitness = float(len(set(scores)) / float(len(self.random_points))) + missing_penalty
		##print(scores)
		#return self.fitness

	#def calc_fitness_oo_na(self):
		#"""
		#Additional penalty if all the partition scores are the same
		#"""
		#scores = []
		#for P in self.random_points:
			#scores.append(self.partition(P))
		#missing_penalty = float((self.L - len(set(scores))) * self.fit_penalty)
		#self.fitness = float(len(set(scores)) / float(self.random_points)) + missing_penalty
		##print(scores)
		#return self.fitness

	#def calc_fitness_mean(self):
		#"""
		#Ol' Fashion Mean
		#"""
		#points = zip(self.init_a,self.init_b)
		#mean_sum = 0
		#for p in points:
			#mean_sum = self.partition(Point(p[0],p[1])) + 1
		#self.fitness = float(mean_sum / self.L)
	
		#return self.fitness

	#def calc_fitness_std(self,weight=1.0):
		#"""
		#Standard deviation of initial guess set, want to maximize
		#"""
		#scores = []
		#for P in self.random_points:
			#scores.append(self.partition(P))
		##print(scores)
		#np.array(scores)
		#return np.std(scores) * weight

	#def calc_fitness_Q_partition(self,weight=1.0):
		#"""
		#Number of partition guess that end in Q partition, want to maximize
		#"""
		#score = 0.0
		#q_part = self.partition(self.Q)
		##print("Q> "+str(q_part))
		#for P in self.random_points:
			##print(self.partition(P))
			#if int(self.partition(P)) == int(q_part):
				#score += 1
		
		#return float(score)/float(self.L) * weight

	#Can't be used with GA rep, it is impossible to determine what Q partition is using pattern
	#def calc_fitness_diff_mean_ga(self,weight=1.0):
		#"""
		#GA Rep, Tries to get mean of all the partition guesses close to partition of Q
		#"""
		#scores = []
		#q_part = self.partition(self.Q)
		##print("Q:"+str(q_part))
		
		#for P in self.random_points:
			##print(self.partition(P), end=" ")
			#scores.append(self.partition(P))
			
		#diff_mean = float(abs(q_part - (sum(scores) / float(len(scores))))/self.L)
		##print("\n"+str(diff_mean))
		#return diff_mean * weight

### GP Fitness ###
### Relating to GP version of the problem
### Using Point Coords and partition(P)
### NTS: will need to eventually postpend _gp to these functions

	def final_fitness_gp(self):
		"""
		GP Representation. Makes a total fitness score
		"""
		#self.fitness = self.calc_fitness_alo()+self.calc_fitness_diff_mean()
		#self.fitness = self.calc_fitness_alo_gp(0.95)+self.calc_fitness_std_gp(0.05)
		#self.fitness = self.calc_fitness_alo_gp()+self.calc_fitness_chisquare_gp(.05)
		#self.fitness = self.calc_fitness_chisquare_gp(0.5)+self.calc_fitness_std_gp(0.5)
		#self.fitness = self.calc_fitness_Q_partition_gp()+self.calc_fitness_match_mean_gp()+self.calc_fitness_penalty_gp()

		#self.fitness = self.calc_fitness_Q_partition_gp(),self.calc_fitness_penalty_gp(thresh=3),
		#self.fitness = self.calc_fitness_match_mean_gp(),self.calc_fitness_Q_partition_gp(),self.calc_fitness_penalty_gp(thresh=2),
		self.fitness = self.calc_fitness_phi_gp(),self.calc_fitness_penalty_gp(thresh=PEN_THRESH),
		return self.fitness,


	def calc_fitness_phi_gp(self,weight=1.0):
		"""
		Using Knuth's result about how approaching Phi makes a hash random
		"""
		
		if FIT_DEBUG: print("Phi_GP")
		scores = []
		phi = (sqrt(5) - 1)/2
		
		#print(self.R)
		
		for r in self.R:
			scores.append(self.partition(r) * phi)
		
		phi_scores = []
		for s in scores:
			phi_scores.append(int(s))
			
		phi_diff = sum(phi_scores) / float(len(phi_scores))
		
		if FIT_DEBUG: print(str(phi_diff))
		return phi_diff * weight


	def calc_fitness_match_mean_gp(self,weight=1.0):
		"""
		GP Representation, based on guess of Q partition, tries to make mean approach that guess
		"""
		if FIT_DEBUG: print("Match_mean_GP")
		scores = []
		for P in self.random_points:
			scores.append(self.partition(P))
			
		match_mean = abs(Q_PARTITION_GUESS - (sum(scores) / float(len(scores))))
		
		if FIT_DEBUG: print(str(match_mean))
		return match_mean * weight


	def calc_fitness_Q_partition_gp(self,weight=1.0):
		"""
		GP Representation, based on guess of Q partition, tries to maximize occurrences of that partition
		"""
		if FIT_DEBUG: print("Q_Part_gp")
		score = 0.0
		q_part = Q_PARTITION_GUESS
		for P in self.random_points:
			if int(self.partition(P)) == int(q_part):
				score += 1
		#q_parts_ratio = 1.0 - float(score)/float(len(self.random_points))
		q_parts_ratio = float(score)/float(len(self.random_points))
		if FIT_DEBUG: print(q_parts_ratio)
		
		return q_parts_ratio * weight

	def calc_fitness_diff_mean_gp(self,weight=1.0):
		"""
		Tries to get mean of all the partition guesses close to partition of Q
		"""
		if FIT_DEBUG: print("Diff_mean_GP")
		scores = []
		for P in self.random_points:
			scores.append(self.partition(P))
			
		diff_mean = float(abs(Q_PARTITION_GUESS - (sum(scores) / float(len(scores)))))
		
		if FIT_DEBUG: print(str(diff_mean))
		return diff_mean * weight


	def calc_fitness_std_gp(self,weight=1.0):
		"""
		GP Representation, Standard Deviation using proper frequency calculation
		"""
		scores = []
		for P in self.random_points:
			scores.append(self.partition(P))

		freqs = []
		for i in range(0,self.L):
			freqs.append(scores.count(i))
		#print(freqs)
		return np.std(freqs) * weight
		
	def calc_fitness_chisquare_gp(self,weight=1.0):
		"""
		Standard Deviation using proper frequency calculation
		"""
		scores = []
		for P in self.random_points:
			scores.append(self.partition(P))

		freqs = []
		for i in range(0,self.L):
			freqs.append(scores.count(i))
		#print(freqs)
		return chisquare(freqs)[0] * weight

	def calc_fitness_alo_gp(self,weight=1.0):
		"""
		At least one value in each partition
		"""
		scores = []
		for P in self.random_points:
			scores.append(self.partition(P))
		score = float(len(set(scores))) / float(len(self.random_points))
		#print(scores)
		return score * weight



	def calc_fitness_penalty_gp(self,thresh=False):
		"""
		GP Rep, Make sure prereqs are there
		"""
		penalty = 0.0	#base penalty
		multi = 1.0	#negative if fitness maximzes
		
		#No coordinate penalty
		#if not re.search('P\.x',self.partition_fn) and not re.search('P\.y',self.partition_fn):
		#	penalty += 100.0

		#No P.x penalty
		#if not re.search('P\.x',self.partition_fn):
		#	penalty += 100.0
		
		#No P.y penalty
		if not re.search('P\.y',self.partition_fn):
			penalty += 10.0

		#Threshold penalty
		scores = []
		for P in self.random_points:
			scores.append(self.partition(P))
		
		if int(thresh):
			if len(set(scores)) <= thresh:
				penalty += 10.0

		#No ints penalty
		#if not re.search('\d',self.partition_fn):
		#	penalty += 100.0

		return multi * penalty

	def partition(self,P):
		"""
		GP Rep, Mutable version of partition fn, always does % self.L
		"""
		#The partition_fn will have values for P.x and P.y which makes this work
		#Best way to handle Points at the Origin?
		# try/catch block is to make sure partition function always works even if mod 0
		
		if P == "Origin":
			#print("origin in partition")
			score = 0
		else:
			try:
				score = int(eval(self.partition_fn))
			except ZeroDivisionError:
				#print("div by zero in partitions ",self.partition_fn)
				score = 1
		return score % self.L

	def fix_fn_notation(self):
		"""
		GP Rep, Changes partition function as generated by DEAP into something 'eval'able
		"""
		self.partition_fn = str.replace(self.partition_fn, 'mul', 'operator.mul')
		self.partition_fn = str.replace(self.partition_fn, 'add', 'operator.add')
		self.partition_fn = str.replace(self.partition_fn, 'neg', 'operator.neg')
		self.partition_fn = str.replace(self.partition_fn, 'sub', 'operator.sub')
		self.partition_fn = str.replace(self.partition_fn, 'protectedMod','operator.mod')
		return True

	def partition_sequence_gp(self):
		"""
		GP Rep, generates a list of current partition function on random points
		"""
		ppath = []
		
		if self.partition_fn == "":
			return False
		
		if self.random_points == []:
			return False
			
		for P in self.random_points:
			ppath.append(self.partition(P))
		
		return ppath

### GA Fit Graveyard
	#def calc_fitness_chisquare_ga(self,weight=1.0):
		#"""
		#Do Chi-square to determine randomness
		#"""
		#scores = []
		#for i in range(0,self.L):
			#scores.append(self.partition_pattern.count(i))
			
		#return chisquare(scores)[0] * weight

	#def calc_fitness_oo_ga(self,weight=1.0):
		#"""
		#Only want 1 per each partition, fit_penalty applied to each non-unique value
		#"""
		#score = 0
		##scores = []
		##for P in self.random_points:
			##scores.append(self.partition(P))
		##missing_penalty = float((self.L - len(set(scores))) * self.fit_penalty)
		##self.fitness = float(len(set(scores)) / float(len(self.random_points))) + missing_penalty
		###print(scores)
		#return score * weight
		


	#def calc_fitness_std_ga(self,weight=1.0):
		#"""
		#GA Rep, Standard deviation of initial guess set, want to maximize
		#"""
		#scores = []
		#for i in range(0,self.L):
			#scores.append(self.partition_pattern.count(i))
		#np.array(scores)
		#return np.std(scores) * weight

### GA Fitness ###
### Relating to GA Represenation of the problem
### No Point Coords and partition_next()
### Different Fitness functions used due to different amount of info available

	def final_fitness_ga(self):
		"""
		GA Representation. Makes a total fitness score
		"""
		#self.fitness = self.calc_fitness_alo_ga()+self.calc_fitness_std_ga()
		#self.fitness = self.calc_fitness_std_ga()
		#self.fitness = self.calc_fitness_alo_ga(0.7)+self.calc_fitness_std_ga(0.3)+self.calc_fitness_penalty_ga(thresh=2)
		#self.fitness = self.calc_fitness_match_mean_ga()+self.calc_fitness_penalty_ga(thresh=2)
		#self.fitness = self.calc_fitness_Q_partition_ga()+self.calc_fitness_penalty_ga(thresh=3)
		#self.fitness = self.calc_fitness_std_ga()+self.calc_fitness_Q_partition_ga()
		#self.fitness = self.calc_fitness_alo_ga(),self.calc_fitness_Q_partition_ga()
		
		#self.fitness = self.calc_fitness_Q_partition_ga(),self.calc_fitness_penalty_ga(thresh=3),
		#self.fitness = self.calc_fitness_std_ga(),self.calc_fitness_Q_partition_ga(),self.calc_fitness_alo_ga(),
		#self.fitness = self.calc_fitness_match_mean_ga(),self.calc_fitness_Q_partition_ga(),self.calc_fitness_penalty_ga(thresh=2),
		
		self.fitness = self.calc_fitness_phi_ga(),self.calc_fitness_penalty_ga(thresh=PEN_THRESH),
		return self.fitness,

	#def calc_fitness_mean_ga(self,weight=1.0):
	#	if FIT_DEBUG: print("Mean_GA")
		

	def calc_fitness_phi_ga(self,weight=1.0):
		"""
		Using Knuth's result about how approaching Phi makes a hash random
		"""
		if FIT_DEBUG: print("Phi_GA")
		scores = []
		phi = (sqrt(5) - 1)/2
		
		#for P in self.random_points:
		#	scores.append(self.partition(P) * phi)
		
		for S in self.partition_pattern:
			scores.append(self.partition_next(0) * phi)
		
		
		phi_scores = []
		for s in scores:
			phi_scores.append(int(s))
			
		phi_diff = sum(phi_scores) / float(len(phi_scores))
		
			
		#phi_diff = abs(phi - (sum(scores) / float(len(scores))))
		#phi_diff = (sum(scores) / float(len(scores)))
		
		if FIT_DEBUG: print(str(phi_diff))
		return phi_diff * weight

	def calc_fitness_alo_ga(self,weight=1.0):
		"""
		GA Rep, At least one value in each partition
		"""
		scores = []
		for P in self.partition_pattern:
			scores.append(self.partition_next(P))
		score = float(len(set(scores))) / float(len(self.partition_pattern))
		return score * weight

	def calc_fitness_match_mean_ga(self,weight=1.0):
		"""
		GA Representation, based on guess of Q partition, tries to make mean approach that guess
		"""
		if FIT_DEBUG: print("Match_mean_GA")
		scores = []
		scores = self.partition_pattern
		match_mean = abs(Q_PARTITION_GUESS - (sum(scores) / float(len(scores))))
		
		if FIT_DEBUG: print(str(match_mean))
		return match_mean * weight

	def calc_fitness_Q_partition_ga(self,weight=1.0):
		"""
		GA Representation, based on guess of Q partition, tries to maximize occurrences of that partition
		expressed as 1-ratio to get minimizable value
		"""
		if FIT_DEBUG: print("Q_Part_GA")
		score = 0.0
		q_part = Q_PARTITION_GUESS
		for P in self.partition_pattern:
			if int(self.partition_next(P)) == int(q_part):
				score += 1
		q_parts_ratio = float(score)/float(len(self.partition_pattern))
		if FIT_DEBUG: print(q_parts_ratio)
		return q_parts_ratio * weight
		
	def calc_fitness_std_ga(self,weight=1.0):
		"""
		GA Representation, Standard Deviation using proper frequency calculation
		"""
		if FIT_DEBUG: print("STD_GA")
		scores = self.partition_pattern

		freqs = []
		for i in range(0,self.L):
			freqs.append(scores.count(i))
		if FIT_DEBUG: print(freqs)
		
		std = np.std(freqs)
		if FIT_DEBUG: print(std)
		return std * weight

	def calc_fitness_penalty_ga(self,weight=1.0,thresh=False):
		"""
		GA Rep, Make sure prereqs are there
		"""
		penalty = 0.0	#base penalty
		multi = 1.0	#positive if fitness maximzes
		
		if FIT_DEBUG: print("Penalty_Ga")
		
		#Clumped values penalty
		scores = self.partition_pattern
		if int(thresh):
			if len(set(scores)) <= thresh:
				penalty += 10.0
		
		if FIT_DEBUG: print(penalty)
		
		return multi * penalty

	def partition_next(self,P):
		"""
		GA Rep, Picks next position in the partition_pattern, 0 indexed
		"""
		#The next parition we pick doesn't rely on the point given to it so we just ignore it
		# pattern_next
		
		#Basic check
		if len(self.partition_pattern) == 0:
			return False
		next_guess = self.partition_pattern[self.pattern_next]
		self.pattern_next = (self.pattern_next + 1) % len(self.partition_pattern)
		
		return int(next_guess)

	def partition_sequence_ga(self):
		"""
		GA Rep generates a list of current partition function on random points
		"""
		ppath = []
		
		if len(self.partition_pattern) == 0:
			return False
		
		#if self.random_points == []:
		#	return False
			
		#this is exactly the partition_pattern, just print that
		for pp in self.partition_pattern:
			print(pp, end = " ")
		
		
		return ppath

### Helper Functions etc for the Algebra of curve math
###
###
	def to_string(self):
		"""
		prints Curve as an algebraic expression
		"""
		return "F("+str(self.prime)+"): y^2 = x^3 + "+str(self.a)+"*x + " + str(self.b)

	def points_equal(self,P,Q):
		"""Determines if two points are the same"""
		if not (self.valid(P) and self.valid(Q)):
			raise ValueError("Invalid inputs")
		else:
			if (P == self.Origin) and (Q == self.Origin):
				return True
			elif (P == self.Origin) and (Q != self.Origin):
				return False
			elif (P != self.Origin) and (Q == self.Origin):
				return False
		if P.x == Q.x and P.y == Q.y:
			return True
		else:
			return False

	def orig_points_equal(self,P,Q):
		"""Determines if two points are the same"""
		if not (self.valid(P) and self.valid(Q)):
			raise ValueError("Invalid inputs")
		else:
			if (P == self.Origin) and (Q == self.Origin):
				return True
			elif (P == self.Origin) and (Q != self.Origin):
				return True
			elif (P != self.Origin) and (Q == self.Origin):
				return True
		if P.x == Q.x and P.y == Q.y:
			return True
		else:
			return False

	def valid(self,P):
		"""
		Calculates if point is actually on Curve, algebraically
		"""
		if P == self.Origin:
			return True
		else:
			return (
				(P.y**2 - (P.x**3 + self.a*P.x + self.b)) % self.prime == 0 and
				0 <= P.x < self.prime and 0 <= P.y < self.prime)

	def inv_mod_p(self,x):
		"""
		Compute an inverse for x modulo p, assuming that x
		is not divisible by p.
		"""
		if x % self.prime == 0:
			raise ZeroDivisionError("Impossible inverse")
		return pow(x, self.prime-2, self.prime)

	def inv(self,P):
		"""
		Inverse of the point P on the elliptic curve y^2 = x^3 + ax + b.
		"""
		if P == self.Origin:
			return P
		return Point(P.x, (-P.y)%self.prime)

	def add(self,P, Q):
		"""
		Sum of the points P and Q on the elliptic curve y^2 = x^3 + ax + b.
		"""
		if not (self.valid(P) and self.valid(Q)):
			raise ValueError("Invalid inputs")

		# Deal with the special cases where either P, Q, or P + Q is
		# the origin.
		if P == self.Origin:
			result = Q
		elif Q == self.Origin:
			result = P
		elif Q == self.inv(P):
			result = self.Origin
		else:
			# Cases not involving the origin.
			if P == Q:
				dydx = (3 * P.x**2 + self.a) * self.inv_mod_p(2 * P.y)
			else:
				dydx = (Q.y - P.y) * self.inv_mod_p(Q.x - P.x)
			x = (dydx**2 - P.x - Q.x) % self.prime
			y = (dydx * (P.x - x) - P.y) % self.prime
			result = Point(x, y)

		# The above computations *should* have given us another point
		# on the curve.
		assert self.valid(result)
		return result

	def mul(self,P,S):
		"""
		Multiplies point S times. ie P + P + P 
		"""
		result = self.add(P,P)
		for s in range(2,S,1):
			result = self.add(result,P)
			
		assert self.valid(result)
		return result

	def rand_point_exhaustive(self):
		"""
		Calculates a random coordinate in the Field for this curve
		"""
		while True:
			x = randint(0,int(self.prime))
			x_side_mod = (x**3 + self.a*x + self.b) % self.prime
			yset = range(0,int(self.prime))
			for y in yset:
				y_side_mod = (y**2) % self.prime
				if x_side_mod == y_side_mod:
					if bool(getrandbits(1)):
						return Point(x, (self.prime - y) % self.prime)
					else:
						return Point(x, (self.prime + y) % self.prime)

	def rand_point(self):
		"""
		Calculates a random point by multiplying a point randomly
		"""
		
		while True:
			mul_scale = randint(0,MAXPOINTRANGE)
			if mul_scale != 1000:
				break
		
		return self.mul(self.P,mul_scale)

	def create_point(self,x,y):
		"""
		Create point on this curve using params, make sure it is valid
		""" 
		p = Point(x,y)
		if self.valid(p):
			return p
		else:
			return False

	def calc_all_points(self):
		"""
		Calculates all of the points in the Field for Given curve
		"""
		points = []
		for x in range(0,int(self.prime)):
			x_side_mod = (x**3 + self.a*x + self.b) % self.prime
			#print x_side_mod
			for y in range(0,int(self.prime)):
				y_side_mod = (y**2) % self.prime
				if x_side_mod == y_side_mod:
					self.points.append(Point(x, (self.prime - y) % self.prime))
					self.points.append(Point(x, (self.prime + y) % self.prime))

		[self.points.append(item) for item in points if item not in self.points]
		self.points = list(set(self.points))
		# +1 for the point at infinity
		self.cardinality = len(self.points) + 1

	def find_k_via_mul(self,P,Q):
		"""Find k*P = Q by adding P + P + P ... until k"""
		k = 1
		temp_point = Point(P.x,P.y)
		while not self.points_equal(temp_point,Q):
			k += 1
			temp_point = self.add(temp_point,P)
		return k

### Generates collection of meaty original runs
###
###
def build_above():
	"""
	Generates collection of meaty original runs
	"""
	o_file = open(str(TC[3])+"_orig.txt","w")
	thresh = int(1.0308 * sqrt(TC[3]))
	print(thresh)
	count = 0
	seed = 3
	while(count <= 1):
		score = rho_set_seed(seed)
		seed += 1
		print(str(count)+","+str(seed)+","+str(score))
		if (score >= thresh):
			score = 0
			print(score)
			o_file.write(str(score)+","+str(seed)+"\n")
			o_file.flush()
			count += 1
		
	o_file.close()
	
### ECCurve object is closed
### These are helper functions to prep/test curve info
###
def generate_long_runs():
	"""
	Will randomly generate rho runs and count of many iterations they take
	"""
	for run in range(0,60):
		seed = randint(0,5000)
		i,s = rho_set_seed(seed)
		longrunf = open(LONGRUN_FILE,"a")
		longrunf.write(i,s)
		longrunf.close()
		
	return True
		
def rho_set_seed(seed):
	"""
	Used for either settling on a seed or to generate multiple original runs
	"""
	random.seed(seed)
	#print("Exploring..."+str(TC[0])+", s "+str(seed)+", l "+str(PARTITIONS)+",",end= " ")

	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.L = PARTITIONS
	prime_order = TC[3]
	
	eC.init_c_prime = INIT_C
	eC.init_d_prime = INIT_D


	for a in range(0,PARTITIONS):
		eC.init_a.append(randint(0,int(TC[3])))

	for a in range(0,PARTITIONS):
		eC.init_b.append(randint(0,int(TC[3])))
		
	eC.P = eC.create_point(TC[4],TC[5])
	eC.Q = eC.create_point(TC[6],TC[7])
	
	#eC.build_R()
	eC.partition_fn = "(P.x)"
	
	k = eC.rho(prime_order,True)
	
	#print(eC.init_a)
	#print(eC.init_b)
	
	print(str(eC.rho_iters)+","+str(seed))
	return(str(eC.rho_iters),str(seed))
	
	
def rho_gp_with_params():
	"""
	Runs a to completion GP RHO using what is in curve_parameters.py
	no cache checking
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.L = PARTITIONS
	prime_order = TC[3]
	
	eC.init_a = MUL_A
	eC.init_b = MUL_B
	eC.init_c_prime = INIT_C
	eC.init_d_prime = INIT_D
	
	eC.P = eC.create_point(TC[4],TC[5])
	eC.Q = eC.create_point(TC[6],TC[7])
	
	eC.partition_fn = FORCE_PARTITION_FUNCTION
	k = eC.rho(prime_order,True)
	print(str(eC.rho_iters))
	return True
	
def rho_ga_with_params():
	"""
	Runs a to completion GA RHO using what is in curve_parameters.py
	no cache checking
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.L = PARTITIONS
	prime_order = TC[3]
	
	eC.init_a = MUL_A
	eC.init_b = MUL_B
	eC.init_c_prime = INIT_C
	eC.init_d_prime = INIT_D
	
	eC.P = eC.create_point(TC[4],TC[5])
	eC.Q = eC.create_point(TC[6],TC[7])
	
	eC.partition_pattern = FORCE_PARTITION_PATTERN
	k = eC.rho_ga(prime_order,True)
	print(str(eC.rho_iters))
	return True

def rho_with_uber_seed():
	"""
	Perform Original Rho with UBER_SEED set int curve_parameters.py
	"""
	random.seed(UBER_SEED)

	if CURVE_DEBUG: print(str(TC[0])+", s "+str(UBER_SEED)+", l "+str(PARTITIONS)+",",end= " ")

	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.L = PARTITIONS
	prime_order = TC[3]
	
	eC.init_a = MUL_A
	eC.init_b = MUL_B
	eC.init_c_prime = INIT_C
	eC.init_d_prime = INIT_D
	
	eC.P = eC.create_point(TC[4],TC[5])
	eC.Q = eC.create_point(TC[6],TC[7])
	
	eC.partition_fn = FORCE_PARTITION_FUNCTION
	
	k = eC.rho(prime_order,True)
	
	ofile = open(ORIGINAL_FILE+"_s_"+str(UBER_SEED)+"_l"+str(eC.L)+".txt","w")
	ofile.write("Original Pollard Rho for: "+str(TC[0])+"\n\n")
	ofile.write(eC.to_string()+"\n")
	ofile.write(str(eC.partition_fn)+"\n")
	ofile.write("Seed: "+str(UBER_SEED)+"\n")
	ofile.write("L "+str(eC.L)+"\n")
	ofile.write("P "+str(eC.P)+"\n")
	ofile.write("Q "+str(eC.Q)+"\n")
	ofile.write("Expected Runs:   "+str(int(sqrt(prime_order)))+"\n")
	ofile.write("Actual Runs:     "+str(eC.rho_iters)+"\n")
	
	if CURVE_DEBUG: print(" > "+str(eC.rho_iters))
	return True

def brute_force_ga():
	"""
	Will do every permutation of partition pattern to determine best results
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.L = PARTITIONS
	prime_order = TC[3]
	
	eC.init_a = MUL_A
	eC.init_b = MUL_B
	eC.init_c_prime = INIT_C
	eC.init_d_prime = INIT_D

	eC.P = eC.create_point(TC[4],TC[5])
	eC.Q = eC.create_point(TC[6],TC[7])
	
	
	for pp in itertools.combinations_with_replacement([7,6,5,4,3,2,1,0],PARTITION_PATTERN_LEN):
		eC.partition_pattern = list(pp)
		eC.rho_ga(prime_order,True)
		print(eC.rho_iters,",",list(pp))

##
# GP Cache Functions
def check_evo_cache_gp(candidate_sol):
	"""
	Checks to see if given pattern is in cache, returns rho
	iters if found
	"""
	try:
		cf = open(GP_EVO_CACHE_BASE,"r")
		for c in cf.readlines():
			match_pattern = list(map(int,c.split(",")[1:-1]))
			#print(match_pattern)
			if match_pattern == candidate_sol:
				print("Match")
				return int(c.split(",")[0])
		return False
	except:
		cf = open(GP_EVO_CACHE_BASE,"a")
		return False

def add_evo_cache_gp(rho_iters, candidate_sol):
	"""
	Add this newly found candidate solution to the cache
	"""
	c_st = ""
	cf = open(GP_EVO_CACHE_BASE,"a")
	for c in candidate_sol:
		c_st += str(c)+","
	cf.write(str(rho_iters)+","+c_st+"\n")
	cf.close()
	return True

def print_evo_cache_gp():
	"""
	Prints contents of cache to screen
	"""
	cf = open(GP_EVO_CACHE_BASE,"r")
	for c in cf.readlines():
		print(c)
	cf.close()
	return True

##
# GA Cache Functions
def check_evo_cache_ga(candidate_sol):
	"""
	Checks to see if given pattern is in cache, returns rho
	iters if found
	"""
	try:
		cf = open(GA_EVO_CACHE_BASE,"r")
		for c in cf.readlines():
			match_pattern = list(map(int,c.split(",")[1:-1]))
			#print(match_pattern)
			if match_pattern == candidate_sol:
				print("Match")
				return int(c.split(",")[0])
		return False
	except:
		cf = open(GA_EVO_CACHE_BASE,"a")
		return False

def add_evo_cache_ga(rho_iters, candidate_sol):
	"""
	Add this newly found candidate solution to the cache
	"""
	c_st = ""
	cf = open(GA_EVO_CACHE_BASE,"a")
	for c in candidate_sol:
		c_st += str(c)+","
	cf.write(str(rho_iters)+","+c_st+"\n")
	cf.close()
	return True

def print_evo_cache_ga():
	"""
	Prints contents of cache to screen
	"""
	cf = open(GA_EVO_CACHE_BASE,"r")
	for c in cf.readlines():
		print(c)
	cf.close()
	return True

## 
# Fitness Snapshots
#
def fitness_snapshot_gp():
	"""
	Shows fitness of initial random points for GP representation
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.L = PARTITIONS
	prime_order = TC[3]
	
	#eC.init_a = MUL_A
	#eC.init_b = MUL_B
	
	#Hopefully this is the same
	for a in range(0,PARTITIONS):
		eC.init_a.append(randint(0,int(TC[3])))

	for a in range(0,PARTITIONS):
		eC.init_b.append(randint(0,int(TC[3])))
	
	#for j in range(0,eC.L):
	#	eC.R.append(eC.add(eC.mul(eC.P,eC.init_a[j]),eC.mul(eC.Q,eC.init_b[j])))

	eC.P = eC.create_point(TC[4],TC[5])
	eC.Q = eC.create_point(TC[6],TC[7])

	eC.init_c_prime = INIT_C
	eC.init_d_prime = INIT_D
	
	eC.partition_fn = FORCE_PARTITION_FUNCTION

	eC.build_R()

	for x,y in zip(TEST_POINTS_X,TEST_POINTS_Y):
		eC.random_points.append(Point(x,y))


	print("GP Rep\n"+eC.to_string())
	print(eC.partition_fn)
	print(eC.final_fitness_gp())
	print(eC.partition_sequence_gp())

def fitness_snapshot_ga():
	"""
	Shows fitness of initial random points for GA representation
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.L = PARTITIONS
	prime_order = TC[3]
	
	eC.init_a = MUL_A
	eC.init_b = MUL_B
	eC.init_c_prime = INIT_C
	eC.init_d_prime = INIT_D
	
	eC.P = eC.create_point(TC[4],TC[5])
	eC.Q = eC.create_point(TC[6],TC[7])
	
	eC.partition_pattern = FORCE_PARTITION_PATTERN

	#Probably don't need this... but negligable really
	#for x,y in zip(TEST_POINTS_X,TEST_POINTS_Y):
	#	eC.random_points.append(Point(x,y))


	print("GA Rep\n"+eC.to_string())
	print(eC.partition_pattern)
	print(eC.final_fitness_ga())

	
def calculate_Q():
	"""
	When we need to find 1000*P
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	P = eC.create_point(TC[4],TC[5])
	Q = eC.mul(P,1000)

	print(str(P.x)+","+str(P.y)+","+str(Q.x)+","+str(Q.y))
	return

def calculate_Q_partition():
	"""
	When we need to find Q partition
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.Q = eC.create_point(int(TC[6]),int(TC[7]))
	eC.L = PARTITIONS
	eC.partition_fn = "(P.x)"
	print(eC.partition(eC.Q))

def random_points_to_file():
	"""
	To speed up process a tad precompute the random set of points
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.P =  eC.create_point(TC[4],TC[5])
	xof = open(POINTS_FILE_X+str(UBER_SEED)+".txt","w")
	yof = open(POINTS_FILE_Y+str(UBER_SEED)+".txt","w")

	for i in range(0,NUM_TEST_POINTS):
		if CURVE_DEBUG: print(str(i)+": ",end="")
		rp = eC.rand_point()
		if CURVE_DEBUG: print(str(rp))
		xof.write(str(rp.x)+"\n")
		yof.write(str(rp.y)+"\n")

		#xof.write("\t"+str(rp.x)+",\n")
		#yof.write("\t"+str(rp.y)+",\n")
		xof.flush()
		yof.flush()
	xof.close()
	yof.close()
	
	return True

def multi_to_file():
	"""
	generates multipliers for rho work
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	af = open(MUL_A_FILE+str(UBER_SEED)+".txt","w")
	bf = open(MUL_B_FILE+str(UBER_SEED)+".txt","w")

	#A
	for a in range(0,PARTITIONS):
		af.write("\t"+str(randint(0,int(TC[3])))+",\n")
	#B
	for b in range(0,PARTITIONS):
		bf.write("\t"+str(randint(0,int(TC[3])))+",\n")

	return True
	
def prime_to_file():
	"""
	Calculates and saves initial C_prime and D_prime as 1/4 and 3/4 into range
	"""
	cf = open(C_PRIME_FILE,"w")
	df = open(D_PRIME_FILE,"w")
	
	cf.write(str(int(TC[3]*0.25)))
	df.write(str(int(TC[3]*0.75)))

	return True

def curve_prep_to_file():
	"""
	generates points,mul data plus what else for curve
	"""
	if CURVE_DEBUG: print("Prepping curve: "+str(TC[0])+" seeded with "+str(UBER_SEED))
	random.seed(UBER_SEED)

	#if CURVE_DEBUG: print("Gen prime values")
	#prime_to_file()

	#if CURVE_DEBUG: print("Gen multipliers")
	#multi_to_file()

	if CURVE_DEBUG: print("Gen points")
	random_points_to_file()
	

def test_ga_rep():
	"""
	Will test various bit and bobs for GA rep of this problem
	"""
	eC = ECCurve(TC[0],TC[1],TC[2])
	eC.L = PARTITIONS
	prime_order = TC[3]
	
	eC.init_a = MUL_A
	eC.init_b = MUL_B
	eC.init_c_prime = INIT_C
	eC.init_d_prime = INIT_D
	
	eC.P = eC.create_point(TC[4],TC[5])
	eC.Q = eC.create_point(TC[6],TC[7])

	#eC.partition_fn = "(P.x)"
	#print("using \n"+str(eC.partition_fn))
	#k = eC.rho(prime_order,True)
	#print("sol: "+str(k))

	eC.partition_pattern = FORCE_PARTITION_PATTERN
	#print("using \n"+str(eC.partition_pattern))
	#k = eC.rho_ga(prime_order,True)
	#print("sol: "+str(k))



if __name__ == "__main__":
	
	#For a really fresh curve need Q values
	#calculate_Q([])

	if len(sys.argv) == 1:
		print(
	"""usage: 
	a - 30 runs above expected
	b - brute force ga
	c - display cache
	f - fitness snapshot GP style
	g - fitness snapshot GA style
	l - generate some long runs
	
	p - write out mutipliers/points based on UBER_SEED
	q - find value of point Q
	r - rho original, via GP to file, using UBER_SEED
	
	s - complete GP RHO using curve_parameters details
	t - complete GA RHO using curve_parameters details

	anything else - going fishing\n""")
	elif sys.argv[1] == 'a':
		build_above()
	elif sys.argv[1] == 'b':
		brute_force_ga()
	elif sys.argv[1] == 'c':
		print_evo_cache_gp()
		print_evo_cache_ga()
	elif sys.argv[1] == 'p':
		curve_prep_to_file()
	elif sys.argv[1] == 'f':
		fitness_snapshot_gp()
	elif sys.argv[1] == 'g':
		fitness_snapshot_ga()
	elif sys.argv[1] == 'r':
		rho_with_uber_seed()
	elif sys.argv[1] == 'q':
		calculate_Q()
	elif sys.argv[1] == 's':
		rho_gp_with_params()
	elif sys.argv[1] == 't':
		rho_ga_with_params()
	elif sys.argv[1] == 'l':
		generate_long_runs()
	else:
		rho_set_seed(sys.argv[1])

	
	#signal_complete()
