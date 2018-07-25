# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:14:43 2016
rotate_morphology
@author: superuser
"""
f=open('morphology_files/T22N5.hoc')
outFile = open('morphology_files/T22N5_try.hoc','a+')
for line in f:
   line = line.split()
   if line[0][0:7]=='pt3dadd': # rotate about y axis for 180 degrees
       outLine = line[0][0:8] + '-1*' + line[0][8::] + line[1] + '-1*' + line[2] + line[3] + '\n'
   else:
       "".join(line)
       outLine = line
   outFile.write(outLine)

outFile.close()
