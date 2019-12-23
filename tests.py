import unittest
from astropy.io import fits
import numpy as np
import os
import grav.grav_fly as grav_fly

class gravTestCase(unittest.TestCase):
	def test_int(self):
		res = grav_fly.main(530000,20,0)

		for i in res:
			self.assertFalse(i is None)


if __name__ == '__main__':
	unittest.main()
