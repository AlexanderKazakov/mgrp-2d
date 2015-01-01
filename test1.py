#!/usr/bin/python
import commands
import os
import shutil
import threading
import sys
import math
import time

# Script for testing mgrp-2d on stupid errors


taskFiles = ["OneStraightDiagonalFractureConstPressure", \
			"OneCurveFracturePolynomialPressure", \
			"ThreeHorizontalFracturesConstPressureSimpleRotation", \
			"ThreeHorizontalFracturesConstPressure", \
			"ThreeVerticalFracturesPolynomialPressure", \
			"ThreeVerticalFracturesConstPressure", \
			"FourFracturesAlongWellConstPressure", \
			"FourFracturesAlongWellPolynomialPressure"]

for taskFile in taskFiles:
	os.system('rm -f results/test1/' + taskFile + '.png')

for taskFile in taskFiles:
	os.system('./mgrp-2d -t ' + 'tasks/tests/' + taskFile + '.xml')
	os.system('cp fractures.png results/test1/' + taskFile + '.png')

