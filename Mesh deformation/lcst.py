"""
Local-CST aerofoil parameterisation method
Author:			Jed Hollom
Co-Authors:		Feng Zhu
				Gabriele Mura
"""

import numpy as np
from scipy.misc import factorial as fact
import scipy.optimize as optimize


class lcst():
	"""Implementation of the local class shape transformation paramererisation method"""

	def __init__(self, n=11, n1=0.5, n2=1.0, samples=200, exp=0, delta_frac=1):
		self.n1 = n1
		self.n2 = n2
		self.n = int(n)
		self.samples = int(samples)
		self.exp = exp
		self.delta_frac = delta_frac

	@staticmethod
	def _sample_function(x):
		"""Sine function ramping from 0 to 1."""
		return ((np.sin((x * np.pi) - (np.pi / 2))) / 2) + 0.5

	# ----------------------------------------------

	def _bernstein(self, x):
		"""Bernstain polinomial basis"""
		r = np.empty((self.n + 1, len(x)), dtype=x.dtype)
		xcontr = np.linspace(0, 1, self.n + 1)

		for i in range(0, self.n + 1):
			k = fact(self.n) / fact(i) / fact(self.n - i)
			local = self.delta_frac * (np.cos(x - xcontr[i])**self.exp)
			r[i] = (k * (1.0 - x)**(self.n - i) * x**i) * local

		return r

	def _shape_function(self, a, x):
		"""Compute the shape functions"""
		a = np.array(a)
		s = self._bernstein(x)
		sx = np.sum(a[:, None] * s, axis=0)
		return sx

	def _coords(self, x, a, zte, sign=1):
		"""Compute the profile coordinates"""
		cx = x**self.n1 * (1. - x)**self.n2
		sx = self._shape_function(a, x)
		return sign * (cx * sx + zte * x)

	# ------------------------------------------------

	def get_profile(self, a, samp_x_u=[], samp_x_l=[]):
		"""Produce coordinates of profile generated using lcst variables"""
		samples = np.linspace(0, 1, self.samples + 1)

		if len(samp_x_u) == 0:
			samp_x_u = self._sample_function(samples[::-1])
		if len(samp_x_l) == 0:
			samp_x_l = self._sample_function(samples)

		samp_x_u = np.array(samp_x_u)
		samp_x_l = np.array(samp_x_l)

		a = np.array(a)
		zteu = 0
		ztel = 0

		# Get coords
		y_u = self._coords(samp_x_u, a[:(len(a) / 2)], zteu)
		y_l = self._coords(samp_x_l, a[(len(a) / 2):], ztel)

		# Convert to list
		samp_x_u = samp_x_u.tolist()
		samp_x_l = samp_x_l.tolist()
		y_u = y_u.tolist()
		y_l = y_l.tolist()

		# Output type: [np.ndarrays]
		return [samp_x_u, y_u], [samp_x_l, y_l]

	# ----------------------------------------------

	def _error_func(self, *parameters):
		"""Error calculation for least squared fitting"""
		CP = parameters[0]
		x = np.array(parameters[1])
		y = np.array(parameters[2])
		echo = parameters[3]
		a = CP[0:self.n + 1]
		err = y - self._coords(x, a, y[np.argmax(x)])
		if echo is True:
			print '		Fitting error: ', np.sqrt(np.sum(err**2.0))
		return err

	def fit_profile(self, data_u, data_l, echo=False):
		"""Least squared profile fitting for lcst coefficients"""
		# Split data
		xu = data_u[0]
		yu = data_u[1]
		xl = data_l[0]
		yl = data_l[1]
		# Generate zeros
		x0 = np.zeros(self.n + 1)
		#  upper surface
		if echo is True:
			print '	[Upper fitting aerofoil]'
		CP = optimize.leastsq(self._error_func, x0=x0, xtol=1e-10, args=(xu, yu, echo))[0]
		au = CP[0:self.n + 1]
		# Lower surface
		if echo is True:
			print '	[Lower fitting aerofoil]'
		CP = optimize.leastsq(self._error_func, x0=x0, xtol=1e-10, args=(xl, yl, echo))[0]
		al = CP[0:self.n + 1]
		a = np.concatenate((au, al))
		return a

	# ----------------------------------------------

	def zero_der(self, x, a=None):
		"""Compute zeroth derivative"""
		r = np.empty((self.n + 1))
		xcontr = np.linspace(0, 1, self.n + 1)

		if a is None:
			a = [1] * (self.n + 1)

		for i in range(0, self.n + 1):
			k = fact(self.n) / fact(i) / fact(self.n - i)
			r[i] = a[i] * k * ((1.0 - x)**(self.n - i + self.n2)) * (x**(i + self.n1)) * self.delta_frac * ((np.cos(x - xcontr[i]))**self.exp)

		return r

	def first_der(self, x, a=None):
		"""Compute first derivative"""
		r = np.empty((self.n + 1))
		xcontr = np.linspace(0, 1, self.n + 1)

		if a is None:
			a = [1] * (self.n + 1)

		for i in range(0, self.n + 1):
			k = fact(self.n) / fact(i) / fact(self.n - i)

			# 1
			alpha = ((1.0 - x)**(self.n - i + self.n2)) * (x**(i + self.n1))
			alpha_p_1 = (i + self.n1) * (x**(i + self.n1 - 1.0)) * ((1.0 - x)**(self.n - i + self.n2))
			alpha_p_2 = (self.n - i + self.n2) * ((1.0 - x)**(self.n - i + self.n2 - 1.0)) * (x**(i + self.n1))
			alpha_p = alpha_p_1 - alpha_p_2

			# 2
			beta = np.cos(x - xcontr[i])**self.exp
			beta_p = - self.exp * (np.cos(x - xcontr[i])**(self.exp - 1.0)) * np.sin(x - xcontr[i])

			r[i] = a[i] * k * self.delta_frac * ((alpha * beta_p) + (beta * alpha_p))

		return r

	def second_der(self, x, a=None):
		"""Compute second derivative"""
		r = np.empty((self.n + 1))
		xcontr = np.linspace(0, 1, self.n + 1)

		if a is None:
			a = [1] * (self.n + 1)

		for i in range(0, self.n + 1):

			k = fact(self.n) / fact(i) / fact(self.n - i)
			cos = np.cos(x - xcontr[i])
			sin = np.sin(x - xcontr[i])

			# 1 ------------------
			theta_1 = (x**(i + self.n1)) * ((1.0 - x)**(self.n - i + self.n2))
			theta_1_p_1 = (x**(i + self.n1)) * (self.n - i + self.n2) * ((1.0 - x)**(self.n - i + self.n2 - 1.0)) * - 1.0
			theta_1_p_2 = (i + self.n1) * (x**(i + self.n1 - 1.0)) * ((1.0 - x)**(self.n - i + self.n2))
			theta_1_p = theta_1_p_1 + theta_1_p_2

			omega_1 = (cos**(self.exp - 1.0)) * sin
			omega_1_p_1 = (cos**(self.exp - 1.0)) * cos
			omega_1_p_2 = (self.exp - 1.0) * (cos**(self.exp - 1.0)) * sin * sin * -1.0
			omega_1_p = omega_1_p_1 + omega_1_p_2

			# 2 ------------------
			theta_2 = (x**(i + self.n1)) * ((1.0 - x)**(self.n - i + self.n2 - 1.0))
			theta_2_p_1 = (x**(i + self.n1)) * (self.n - i + self.n2 - 1.0) * ((1.0 - x)**(self.n - i + self.n2 - 2.0)) * - 1.0
			theta_2_p_2 = (i + self.n1) * (x**(i + self.n1 - 1.0)) * ((1.0 - x)**(self.n - i + self.n2 - 1.0))
			theta_2_p = theta_2_p_1 + theta_2_p_2

			omega_2 = cos**self.exp
			omega_2_p = - self.exp * (cos**(self.exp - 1.0)) * sin

			# 3 ------------------
			theta_3 = (x**(i + self.n1 - 1.0)) * ((1.0 - x)**(self.n - i + self.n2))
			theta_3_p_1 = (x**(i + self.n1 - 1.0)) * (self.n - i + self.n2) * ((1.0 - x)**(self.n - i + self.n2 - 1.0)) * - 1.0
			theta_3_p_2 = (i + self.n1 - 1.0) * (x**(i + self.n1 - 2)) * ((1.0 - x)**(self.n - i + self.n2))
			theta_3_p = theta_3_p_1 + theta_3_p_2

			omega_3 = cos**self.exp  # SAME AS OMEGA_2
			omega_3_p = - self.exp * (cos**(self.exp - 1.0)) * sin  # SAME AS OMEGA_2_P

			# Summation ----------
			term_1 = self.exp * ((theta_1 * omega_1_p) + (theta_1_p * omega_1))
			term_2 = (self.n - i + self.n2) * ((theta_2 * omega_2_p) + (theta_2_p * omega_2))
			term_3 = (i + self.n1) * ((theta_3 * omega_3_p) + (theta_3_p * omega_3))
			r[i] = a[i] * k * self.delta_frac * (- term_1 - term_2 + term_3)

		return r

	#  ----------------------------------------------

	def get_bernsteins(self, cst_ceofs):
		"""Generate plotting data for bernstein polynomials"""
		x = np.array([(1.0 / self.samples) * item for item in range(0, self.samples + 1)])

		# Get shape function
		dir_0_u = []
		dir_0_l = []
		for entry in x:
			dir_0_u.append(self.zero_der(entry))
			dir_0_l.append(self.zero_der(entry))

		# Scale polynomials
		y_u = []
		y_l = []
		for p in range(0, len(dir_0_u[0])):
			n_u = [i[p] * cst_ceofs[p] for i in dir_0_u]
			n_l = [i[p] * cst_ceofs[(len(cst_ceofs) / 2) + p] for i in dir_0_l]
			y_u.append(n_u)
			y_l.append(n_l)

		return x, y_u, y_l
