#!/usr/bin/python
import commands
import os
import shutil
import threading
import sys
import math
import time

# Script for testing mhf-2d on stupid errors


taskFiles = ["OneStraightDiagonalFractureConstPressure", \
             "OneCurveFracturePolynomialPressure", \
             "TwoFractures", \
             "ThreeHorizontalFracturesConstPressureSimpleRotation", \
             "ThreeHorizontalFracturesConstPressure", \
             "ThreeVerticalFracturesPolynomialPressure"]

for taskFile in taskFiles:
	os.system('rm -f results/FastTest/' + taskFile + '.png')

for taskFile in taskFiles:
	os.system('./mhf-2d -t ' + 'tasks/tests/FastTest/' + taskFile + '.xml')
	os.system('cp fractures.png results/FastTest/' + taskFile + '.png')

