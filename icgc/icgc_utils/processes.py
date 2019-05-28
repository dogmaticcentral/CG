#! /usr/bin/python3
#
# This source code is part of icgc, an ICGC processing pipeline.
# 
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#
import multiprocessing
import os

########################################
def get_process_id():

	return os.getpid()

# don't know how to do this if there are no other_args,
# except by passing it an empty list in th place
########################################
def partition_load(number_of_chunks, load_length):
	partition = {}
	load = []
	for thr in range(number_of_chunks):
		load.append(0)

	for job in range (load_length):
		load[(job%number_of_chunks)] += 1

	total = 0
	for ps in range (number_of_chunks):
		ps_from = total
		ps_to   = total + load[ps]
		total  += load[ps]

		if (ps_from >= load_length):
			break
		if (ps == number_of_chunks-1):
			ps_to = load_length
		partition[ps] = [ps_from, ps_to]
	return partition

########################################
def linear_pll(number_of_chunks, embarassingly_pllbl_fn, input_list, other_args, return_dict=None):

	partition = partition_load(number_of_chunks, len(input_list))

	processes = []
	for ps in range (number_of_chunks):
		[ps_from, ps_to] = partition[ps]
		args = (input_list[ps_from:ps_to], other_args)
		# addint to a tuple
		#without the comma,  () interpreted as the order of precedence brackets
		if return_dict!=None: args+=(return_dict,)
		process = multiprocessing.Process(target=embarassingly_pllbl_fn, args=args)
		try:
			process.start()
			processes.append(process)
		except:
			print("Error: unable to start process")
			return False

	return processes


########################################
def weighted_partition(number_of_chunks, input_list,  weights):
	if len(input_list)!=len(weights):
		print("input list and the weights must be of the same length")
		exit()
	if len(input_list)==0:
		return [],[]
	input_sorted   = sorted(input_list, key=lambda t: weights[input_list.index(t)], reverse=True)
	weights_sorted = sorted(weights, reverse=True)
	input_list = input_sorted
	weights = weights_sorted
	#partition = [[]]*number_of_chunks # bizarre: this is n ptrs to the same list, set to be empty
	partition = []
	for i in range(number_of_chunks): partition.append([])
	partition_weight = [0]*number_of_chunks
	for i, element in enumerate(input_list):
		min_part_wt, min_part_idx = min((val, idx) for (idx, val) in enumerate(partition_weight))
		partition[min_part_idx].append(element)
		partition_weight[min_part_idx] += weights[i]
	return partition, partition_weight

##############
def weighted_pll(number_of_chunks, embarassingly_pllbl_fn, input_list,  weights, other_args, return_dict):
	partition, partition_weight = weighted_partition(number_of_chunks, input_list,  weights)
	# for s, sublist in enumerate(partition):
	# 	print(" ========= sublist  wt: %10d  (%.2e )=======" % (partition_weight[s],partition_weight[s]) )
	# 	for idx, element in enumerate(sublist):
	# 		print(element, weights[idx])
	# exit()
	processes = []
	for ps in range (number_of_chunks):
		args = (partition[ps], other_args)
		if return_dict!=None: args+=(return_dict,)
		process = multiprocessing.Process(target=embarassingly_pllbl_fn, args=args)
		try:
			process.start()
			processes.append(process)
		except:
			print("Error: unable to start process")
			return False

	return processes



########################################
def round_robin_pll(number_of_chunks,embarassingly_pllbl_fn, list, other_args):
	list_per_process = []
	for process in range(number_of_chunks):
		list_per_process.append([])

	for element in list:
		idx = list.index(element)
		list_per_process[idx%number_of_chunks].append(element)

	# run
	processes = []
	for ps in range (number_of_chunks):

		process = multiprocessing.Process(target=embarassingly_pllbl_fn, args=(list_per_process[ps], other_args))
		try:
			process.start()
			processes.append(process)
		except:
			print("Error: unable to start process")
			return False

	return processes


###########
def parallelize(number_of_chunks, embarassingly_pllbl_fn, list, other_args, round_robin=False):

	if number_of_chunks < 1:
		print("number of processes is expected to be >= 1")
		return False

	if number_of_chunks == 1:
		if other_args==None:
			ret = embarassingly_pllbl_fn(list)
		else:
			ret = embarassingly_pllbl_fn(list, other_args)
		return ret

	if round_robin:
		return round_robin_pll(number_of_chunks,embarassingly_pllbl_fn, list, other_args)
	else:
		return linear_pll(number_of_chunks,embarassingly_pllbl_fn, list, other_args)


###########
def pll_w_return(number_of_chunks, embarassingly_pllbl_fn, input_list, other_args, weights=None):
	if number_of_chunks < 1:
		print("number of processes is expected to be >= 1")
		return False

	if number_of_chunks == 1:
		return_dict = {}
		embarassingly_pllbl_fn(input_list, other_args, return_dict)
		return return_dict

	manager = multiprocessing.Manager()
	return_dict = manager.dict()

	if not weights:
		processes = linear_pll(number_of_chunks, embarassingly_pllbl_fn, input_list, other_args, return_dict)
	else:
		processes = weighted_pll(number_of_chunks, embarassingly_pllbl_fn, input_list,  weights, other_args, return_dict)

	wait_join(processes)

	return return_dict

#########################################
def wait_join(processes):
	for process in processes:
		process.join()


def test_fn(list, other_args, return_dict):
	pid = get_process_id()
	new_list = [l+pid for l in list]
	return_dict[pid] = new_list
	return

def main():
	inlist = [-3, -2, -1, 0, 0, 0, 1, 2, 3 ]
	number_of_chunks = 3
	other_args = []
	return_dict = pll_w_return(number_of_chunks, test_fn, inlist, other_args)
	for k, v  in return_dict.items():
		print (k,v)

	print("==============")
	weighted_partition(3, [1,2,3,4,5,6,7,8],  [36, 25, 18, 7, 5, 3, 1, 1])

#########################################
if __name__ == '__main__':
	main()
