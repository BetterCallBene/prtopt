#!/usr/bin/python

import sys
import os
import tempfile
import string
from ConfigParser import SafeConfigParser

import GenerateScript as GS

def init_and_call(ini_name):

	parser = SafeConfigParser()
	apps_path = os.getcwd()
	ini_path = "%s/%s"%(apps_path, ini_name)
	
	parser = SafeConfigParser()
	parser.read(ini_path)
	

	props={ "template_dir": parser.get('template', 'dir'), "template_name" : parser.get('template', 'file'),
			  "src_dir" : parser.get('src', 'dir'), "src_file" : parser.get('src', 'file'), 
			  "maple_src_dir" : parser.get('maple', 'dir'), "maple_src_file" : parser.get('maple', 'file'),
			  "generateFiles" : parser.get('src', 'generateFiles')
	}
	generateOpts = props["generateFiles"]
	generateOpts = generateOpts[1:-1].replace(' ', '').split(',')

	pathTmp = tempfile.gettempdir() + '/%s'

	ind = 0
	tmp_files =[]
	for name in generateOpts:
		tmp_files.append(pathTmp % name)
		ind = ind + 1

	reps = []

	for name, value in parser.items('replacement'):
		reps.append([name, value])
		#print "%s = %s" %(name, value)

	template_path ="%s/%s/%s" % (apps_path, props["template_dir"], props["template_name"])
	src_path = "%s/%s/%s" %(apps_path, props["src_dir"], props["src_file"])
	maple_path = "%s/%s/%s" %(apps_path, props["maple_src_dir"], props["maple_src_file"])

	GS.generate(parser.get('maple', 'app_path'), template_path, src_path, maple_path, tmp_files, reps)

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "ERROR: Anzahl der Parameter ist falsch %d", len(sys.argv)
	else:
		ini_name = sys.argv[1]
		init_and_call(ini_name)