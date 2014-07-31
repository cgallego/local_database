# -*- coding: utf-8 -*-
"""
Master python script to run each module in sequence, run in debug mode

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

from send2_mydatabase import *

if __name__ == '__main__':    
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__)) 
    print path_rootFolder
       
    # Open filename list
    print sys.argv[1]
    file_ids = open(sys.argv[1],"r")
    file_ids.seek(0)

    print sys.argv[2]
    file_muscleVOI = open(sys.argv[2],"r")
    file_muscleVOI.seek(0)

    for k in range(59):
        line = file_ids.readline()
        line_muscleVOI = file_muscleVOI.readline()
    print line
    
    
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
        SeriesID = fileline[6] # corresponds to dynamic sequence
        finding_side = fileline[7]
        Diagnosis = fileline[11]
        T2SeriesID = fileline[12]
        
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        Send2DB = Send()
        img_folder, rowCase = Send2DB.queryDatabase(cond, StudyID, dateID)
        
        #############################
        ###### 2) Get T2 info
        #############################
        #img_folder = 'Z:/Cristina/MassNonmass/mass/'
        [path_T2Series, T2SeriesID] = Send2DB.processT2(T2SeriesID, img_folder, StudyID, DicomExamNumber)
        
        #############################                  
        ###### 3) Check segmentation accuracy with annotations
        #############################
        [series_path, phases_series, lesion3D] = Send2DB.checkSegment(path_rootFolder, fileline, line, cond, StudyID, DicomExamNumber, SeriesID, Lesionfile, path_T2Series)
        
        #############################
        # 4) Extract Lesion and Muscle Major pectoralies signal                                   
        ############################# 
        [T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR, morphoT2features, textureT2features] = Send2DB.T2_extract(T2SeriesID, path_T2Series, line_muscleVOI, lesion3D)
       
        #############################
        ###### 4) Extract Dynamic features
        #############################
        [dyn_inside, dyn_contour] = Send2DB.extract_dyn(series_path, phases_series, Send2DB.lesion3D)
        
        #############################
        ###### 5) Extract Morphology features
        #############################
        morphofeatures = Send2DB.extract_morph(series_path, phases_series, Send2DB.lesion3D)
        
        #############################        
        ###### 6) Extract Texture features
        #############################
        texturefeatures = Send2DB.extract_text(series_path, phases_series, Send2DB.lesion3D)       

        
        #############################
        ###### 8) End and Send record to DB
        #############################
        Send2DB.addRecordDB_lesion(Lesionfile, Send2DB.casesFrame['id'].iloc[0], DicomExamNumber, dateID, Send2DB.casesFrame.iloc[0], finding_side, Send2DB.dataCase.iloc[0], cond, Diagnosis, 
                           lesion_id, BenignNMaligNAnt,  SeriesID, T2SeriesID)
                           
        Send2DB.addRecordDB_features(lesion_id, dyn_inside, dyn_contour, morphofeatures, texturefeatures)   

        Send2DB.addRecordDB_annot(lesion_id, Send2DB.annot_attrib, Send2DB.eu_dist_mkers, Send2DB.eu_dist_seg)

        Send2DB.addRecordDB_T2(lesion_id, T2SeriesID, Send2DB.dataCase.iloc[0], morphoT2features, textureT2features, T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR)
    
        
        ## continue to next case
        line = file_ids.readline()
        print line
        line_muscleVOI = file_muscleVOI.readline()
        
        
    file_ids.close()
    file_muscleVOI.close()