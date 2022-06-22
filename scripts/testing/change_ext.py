#!/usr/bin/python

# Script to change fastq paired read file extensions on all paried reads in a folder
# Nicola Coyle
# 27/04/2022
# For use with dadaist2 ( produce usable Microbiomeanalyist files first go)

# import modules
import os
import re
import os.path
import shutil
import sys, getopt

# set variables
argv=sys.argv[1:]
directory = ''
current_ext1 = ''
new_ext1 = ''
current_ext2 = ''
new_ext2 = ''
outfolder = ''

# optarg
try:
    opts, args = getopt.getopt(argv,"d:c:n:C:N:o:",["dir=","curext1=","newext1=","curext2=","newext2=","outfolder="])
except getopt.GetoptError:
    print ("""Script to change the file extensions for paired fastq files

  how to use:

  change_ext.py -d <directory> -c <current_ext1> -n <new_ext1> -C <current_ext2> -N <new_ext2> -o <outfolder>""")
    sys.exit(2)
for opt, arg in opts:
    print (opt, arg)
    if opt == '-h':
        print( """Script to change the file extensions for paired fastq files

  how to use:

  change_ext.py -d <directory> -c <current_ext1> -n <new_ext1> -C <current_ext2> -N <new_ext2> -o <outfolder>""")
        sys.exit()
    elif opt in ("-d", "--directory"):
        directory = arg
    elif opt in ("-c", "--current_ext1"):
        current_ext1 = arg
    elif opt in ("-n", "--new_ext1"):
        new_ext1 = arg
    elif opt in ("-C", "--current_ext2"):
        current_ext2 = arg
    elif opt in ("-N", "--new_ext2"):
        new_ext2 = arg
    elif opt in ("-o", "--outfolder"):
        outfolder = arg

# check if required arguments supplied
if not directory: print("Directory not supplied") ; exit()
if not current_ext1: print("Old ext R1 not supplied") ; exit()
if not new_ext1: print("New ext R1 not supplied") ; exit()
if not current_ext2: print("Old ext R2 not supplied") ; exit()
if not new_ext2: print("New ext R2 not supplied") ; exit()

if len(os.listdir(directory)) == 0:
    print("Directory %s empty" %(directory))
    sys.exit()

# check first if a new output directory is specified
# if yes - files will be copied before renaming
# if not - files will be renamed in place
if not outfolder:
    print('File extensions will be replaced in current directory')
    for filename in os.listdir(directory): # loop files
        infilename_p = os.path.join(directory,filename)
        dir_name = os.path.dirname(os.path.abspath(infilename_p))# get absolute path
        infilename = os.path.join(dir_name,filename)
        if current_ext1 in filename: # R1 files
             newname = infilename.replace(current_ext1,new_ext1)
        elif current_ext2 in filename: # R2 files
             newname = infilename.replace(current_ext2,new_ext2)
        else:
             print("Error: file extensions %s and %s not found in file name : %s" %(current_ext1, current_ext2,infilename))
             sys.exit()
        print("newname %s " %(newname))
        if not os.path.isfile(newname):
            output = os.rename(infilename, newname)
            print('file %s changed to %s' %(infilename, newname))
        else:
            print("File already exists")

else:
    print ("Output directory: %s\n Files will be copied to this directory with new extenions" %(outfolder))

    # check if output directory exists and make one if not
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
        print("output directory made: %s" %(outfolder))

    for filename in os.listdir(directory):
        infilename_p = os.path.join(directory,filename)
        dir_name = os.path.dirname(os.path.abspath(infilename_p))# get absolute path
        infilename = os.path.join(dir_name,filename)
        copiedname = os.path.join(outfolder,filename)
        if current_ext1 in filename: # R1 files
             newname = copiedname.replace(current_ext1,new_ext1)
        elif current_ext2 in filename: # R2 files
             newname = copiedname.replace(current_ext2,new_ext2)
        else:
             print("Error: file extensions %s and %s not found in file name : %s" %(current_ext1, current_ext2, copiedname))
             sys.exit()
        if not os.path.isfile(newname):
            shutil.copy(infilename, outfolder) # copy file to output directory
            output = os.rename(copiedname, newname) # change extension
            print('file %s changed to %s' %(infilename, newname))
        else:
            print("File already exists")
