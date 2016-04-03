#!/usr/bin/python

import os
import re
import numpy as np
from shutil import move


#MAPLE_PATH = '/Library/Frameworks/Maple.framework/Versions/Current/bin/maple'
#CONFIG_NAME = 'config.py'

def sub_rep(rep_rep, start_index):

	tmp =rep_rep	
	len_rep = len(rep_rep)
	rep_rep = rep_rep[start_index+1:len_rep]
	
	end_index = rep_rep.find('%')
	rep_rep = rep_rep[0:end_index]
	
	index = rep_rep.find('.')
	start_ind =  int(rep_rep[0:index])
	end_ind = int(rep_rep[index+2:end_index])
	
	return [start_ind, end_ind, start_index + end_index + 2, tmp[start_index:start_index + end_index+2]]

	#start_ind  = int(rep_rep[0:index - 1])
	#print start_ind	
#	end_ind = rep_rep[index + 2:end_index]
#	
			

def generate(MAPLE_PATH, template_path, output, maple_path, tmp_files, reps):
	
	print 'Temporaere loeschen'
	
	for filename in tmp_files:
		if os.path.exists(filename):
			os.remove(filename)

	if os.path.exists(output):
		os.remove(output)

	print 'Lese Template ein'

	hTemplate = open(template_path, 'r')
	try:
		text = hTemplate.read()
	finally:
		hTemplate.close()

	print 'Do maple stuff'
	os.system(MAPLE_PATH + ' < ' + maple_path) 

	print 'Ersetzungen durchfuehren'
	i = 0
	for tmp_file in tmp_files:
		text_rep = u'${%d}$' % i
		htmp_file = open(tmp_file, 'r')

		try:
			text_input = htmp_file.read()
			text = text.replace(text_rep, text_input)
			
		finally:
			htmp_file.close()
		i = i + 1

	for rep in reps:
		text = text.replace(rep[0], rep[1])

# pseudo multiarray

	i = 0

	n = 13
	m = 17

	ind = 0
	while i < n:
		j = 0
		while j < m:
			rep = "pd[%d][%d]" % (i, j)
			to = "pd[%d]" % ind
			text = text.replace(rep, to)
			ind = ind + 1
			j = j + 1
		i = i + 1

	
	l = 17
	i = 0
	ind = 0
	while i < n:
		j = 0
		while j < m:
			k = 0
			while k < l:
				rep = "pdd[%d][%d][%d]" % (i, j, k)
				to = "pdd[%d]" % ind
				text = text.replace(rep, to)
				ind = ind + 1
				k = k + 1
			j = j + 1
		i = i + 1

	i = 0
	n = 2*234 + 234 * 234

	while i < n:
		rep = "DMat[%d][0]" %(i)
		to =  "DMat[%d]" %(i)
		text= text.replace(rep, to)
		i = i + 1


	file_write_to = open(output, 'w')
	file_write_to.write(text)
	file_write_to.close()
	

	for filename in tmp_files:
		if os.path.exists(filename):
			os.remove(filename)


#	template_file = open(template_path, 'r')

#	try:
#		text = template_file.read()
#	finally:
#		template_file.close()

#	for tmp_file in tmp_files:











	
