#!/usr/bin/python
import commands
import os
import shutil
import threading
import sys
import math
import time

# Script for testing mgrp-2d on match with analytical solutions


taskFiles = ["StraightFractureConstPressureForDisplacementChecking", \
			"StraightFractureLagPressureForDisplacementChecking"]


for taskFile in taskFiles:
	os.system('rm -f results/test2/' + taskFile + '.png')

for taskFile in taskFiles:
	os.system('./mgrp-2d -t ' + 'tasks/tests/' + taskFile + '.xml' + ' -d 1')
	os.system('cp displacements.png results/test2/' + taskFile + '.png')

