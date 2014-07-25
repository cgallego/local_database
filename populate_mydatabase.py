# -*- coding: utf-8 -*-
"""
Master python script to run each module in sequence

Arguments:
============
sys.argv[1] = input text file with one case per line in the following format:
cond	StudyNumber	DicomExamNumber	LesionID	StudyDate	SeriesID	BreastSide	PathReportID	PathoBreastSide	BenignNMaligNAnt	Diagnosis	T2SeriesID	annotations		


Created on Tue Jul 15 16:48:07 2014
@ author (C) Cristina Gallego, University of Toronto
----------------------------------------------------------------------
 """
import os, os.path
import sys
from sys import argv, stderr, exit
import numpy as np
import dicom
import psycopg2
import pandas as pd

from send2_mydatabase import *

if __name__ == '__main__':    
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__))
    print path_rootFolder
       
    # Open filename list
    print sys.argv[1]
    file_ids = open(sys.argv[1],"r")
    file_ids.seek(0)
    line = file_ids.readline()
    
    print sys.argv[2]
    file_muscleVOI = open(sys.argv[2],"r")
    file_muscleVOI.seek(0)
    line_muscleVOI = file_muscleVOI.readline()
    
    while ( line ) : 
        # Get the line: StudyNumber    DicomExamNumber    MRN    chosen_lesions_id    StudyDate    SeriesID    image_pos_pat    image_ori_pat
        fileline = line.split()
        lesion_id = int(fileline[0] )
        cond = fileline[1] 
        BenignNMaligNAnt = cond[-1]
        StudyID = fileline[2]  
        DicomExamNumber = fileline[3]
        Lesionfile = fileline[4]
        dateID = fileline[5]
        SeriesID = fileline[6] # corresponds to dynamic sequence;
        finding_side = fileline[7]
        Diagnosis = fileline[11]
        T2SeriesID = fileline[12]
        
        Send2DB = Send()
        Send2DB.process(path_rootFolder, fileline, line, line_muscleVOI, lesion_id, cond, BenignNMaligNAnt, StudyID, DicomExamNumber, Lesionfile, dateID, SeriesID, finding_side, Diagnosis, T2SeriesID)
        
        ## continue to next case
        line = file_ids.readline()
        line_muscleVOI = file_muscleVOI.readline()
        
    file_ids.close()
    file_muscleVOI.close()