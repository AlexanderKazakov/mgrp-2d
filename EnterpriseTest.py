#!/usr/bin/python
import commands
import os
import shutil
import threading
import sys
import math
import time

# Script for enterprise demonstration of mhf-2d


taskFiles = ["ThreeVerticalFracturesConstPressure", \
             "FourFracturesAlongWellConstPressure", \
             "FourFracturesAlongWellPolynomialPressure"]

for taskFile in taskFiles:
	os.system('rm -f results/EnterpriseTest/' + taskFile + '.png')

for taskFile in taskFiles:
	os.system('./mhf-2d -t ' + 'tasks/tests/EnterpriseTest/' + taskFile + '.xml')
	os.system('cp fractures.png results/EnterpriseTest/' + taskFile + '.png')

