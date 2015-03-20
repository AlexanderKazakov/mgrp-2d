#!/usr/bin/python
import commands
import os
import shutil
import threading
import sys
import math
import time

# Script for testing mhf-2d on match with analytical solutions


taskFiles = ["StraightFractureConstPressure135Tip", \
             "StraightFractureLagPressure135Tip", \
             "StraightFractureConstPressureConstTip", \
             "StraightFractureLagPressureConstTip"]


for taskFile in taskFiles:
	os.system('rm -f results/TestOnDisplacementDiscontinuities/' + taskFile + '.png')

for taskFile in taskFiles:
	os.system('./mhf-2d -t ' + 'tasks/tests/TestOnDisplacementDiscontinuities/' + taskFile + '.xml' + ' -d 1')
	os.system('cp displacements.png results/TestOnDisplacementDiscontinuities/' + taskFile + '.png')

