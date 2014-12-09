#!/usr/bin/python
import commands
import os
import shutil
import threading
import sys
import math
import time

# Script for testing mgrp-2d on errors


taskFiles = ["OneFracture", "OneDiagonalFracture", "ThreeHorizontalFractures", "ThreeVerticalFractures"]

for taskFile in taskFiles:
	os.system('./mgrp-2d -t ' + 'tasks/tests/' + taskFile + '.xml')
	os.system('cp fractures.png results/test1/' + taskFile + '.png')

