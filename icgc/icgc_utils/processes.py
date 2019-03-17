
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
def parallelize (number_of_chunks, embarassingly_pllbl_fn, list, other_args):


	if (number_of_chunks < 1):
		print("number of processes is expected to be >= 1")
		return False

	if (number_of_chunks == 1):
		if other_args==None:
			ret = embarassingly_pllbl_fn (list)
		else:
			ret = embarassingly_pllbl_fn (list, other_args)
		return ret


	#########################################
	# nontrivial
	load = []
	for thr in range(number_of_chunks):
		load.append(0)

	for job in range (len(list)):
		load[(job%number_of_chunks)] += 1

	# run
	total = 0
	processes = []
	for ps in range (number_of_chunks):
		ps_from = total
		ps_to   = total + load[ps]
		total  += load[ps]

		if (ps_from >= len(list)):
			break
		if (ps == number_of_chunks-1):
			ps_to = len(list)

		process = multiprocessing.Process(target=embarassingly_pllbl_fn, args=(list[ps_from:ps_to], other_args))
		try:
			process.start()
			processes.append(process)
		except:
			print("Error: unable to start process")
			return False

	return processes


#########################################
def wait_join(processes):
	for process in processes:
		process.join()
