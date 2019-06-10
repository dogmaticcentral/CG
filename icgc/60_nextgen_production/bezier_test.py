#! /usr/bin/python3
from matplotlib import pyplot as plt
import bezier
import numpy as np


#######################
def main():

	# x = [0.0, 1.0, 2.0, 3.0]
	# y = [0.0, 1.0, 2.0, 3.0]
	# ps = [1.0, 2.0]

	# this will crash with
	# File "src/bezier/_speedup.pyx", line 485, in bezier._speedup.curve_intersections
	# ValueError: Curve intersection failed to converge to approximately linear subdivisions after 20 iterations.
	x = [0.0,   69708.0,  204470.0, 339500.0]
	y = [0.000,   0.052,     0.141,  0.214]
	ps = [237685.0, 161745.0]

	# this, however, works
	norm = x[-1]
	x = [e/norm for e in x]
	ps = [p/norm for p in ps]


	nodes = np.asfortranarray([x,y])
	curve = bezier.Curve(nodes, degree=2)

	plt.scatter(x, y)
	plt.title("test")
	plt.ylabel('x')
	plt.xlabel('y')
	ax = plt.gca()
	curve.plot(num_pts=256, ax=ax, color = 'red')
	plt.ylim([0.0, y[-1]*1.1])

	vertical_lines = []
	for p in  ps:
		vertical_lines.append(bezier.Curve(np.asfortranarray([[float(p), float(p)],[0.0, 5.0]]), degree=1))

	for vl in vertical_lines:
		vl.plot(num_pts=25, ax=ax, color = 'red')
	plt.show()

	for vl in vertical_lines:
		intersections = curve.intersect(vl)
		s_vals = np.asfortranarray(intersections[0, :])
		point = curve.evaluate_multi(s_vals)
		if point.size==2:
			print(" >> ", point[1][0])
		else:
			print(" >> ", "no return")






#########################################
if __name__ == '__main__':
	main()
