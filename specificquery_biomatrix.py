# -*- coding: utf-8 -*-
"""
Created to run a specific query on a batch list

Created on Tue Sep 23 13:12:35 2014

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

from sqlalchemy import Column, Integer, String
from sqlalchemy.orm import sessionmaker
import datetime
import database
from base import Base, engine
import pandas as pd


def queryDatabase(session, StudyID, redateID):
    """
    run : Query by StudyID/AccesionN pair study to local folder NO GRAPICAL INTERFACE. default print to output
    
    Inputs
    ======
    StudyID : (int)    CAD StudyID
    redateID : (int)  CAD StudyID Data of exam (format yyyy-mm-dd)
    
    Output
    ======
    """                   
    datainfo = []; 
       
    for cad, exam, in session.query(database.Cad_record, database.Exam_record).\
                filter(database.Cad_record.pt_id==database.Exam_record.pt_id).\
                filter(database.Cad_record.cad_pt_no_txt == str(StudyID)).\
                filter(database.Exam_record.exam_dt_datetime == str(redateID)).all():
                    
        # print results
        if not cad:
            print "cad is empty"
        if not exam:
            print "exam is empty"          
    
        datainfo.append( [([cad.cad_pt_no_txt, cad.latest_mutation_status_int,
                                exam.exam_dt_datetime, exam.mri_cad_status_txt, exam.comment_txt,
                                exam.sty_indicator_rout_screening_obsp_yn,
                                exam.sty_indicator_high_risk_yn,
                                exam.sty_indicator_high_risk_brca_1_yn,  
                                exam.sty_indicator_high_risk_brca_2_yn, 
                                exam.sty_indicator_high_risk_brca_1_or_2_yn, 
                                exam.sty_indicator_high_risk_at_yn, 
                                exam.sty_indicator_high_risk_other_gene_yn, 
                                exam.sty_indicator_high_risk_prior_high_risk_marker_yn, 
                                exam.sty_indicator_high_risk_prior_personal_can_hist_yn, 
                                exam.sty_indicator_high_risk_hist_of_mantle_rad_yn, 
                                exam.sty_indicator_high_risk_fam_hist_yn, 
                                exam.sty_indicator_add_eval_as_folup_yn, 
                                exam.sty_indicator_folup_after_pre_exam_yn,
                                exam.sty_indicator_pre_operative_extent_of_dis_yn,
                                exam.sty_indicator_post_operative_margin_yn, 
                                exam.sty_indicator_pre_neoadj_trtmnt_yn,
                                exam.sty_indicator_prob_solv_diff_img_yn, 
                                exam.sty_indicator_scar_vs_recurr_yn,
                                exam.sty_indicator_folup_recommend_yn, 
                                exam.sty_indicator_prior_2_prophy_mast_yn])] )
        
    ################### Send to table display  
    # add main CAD record table       
    colLabels = ("cad.cad_pt_no_txt", "cad.latest_mutation", "exam.exam_dt_datetime", "exam.mri_cad_status_txt", "exam.comment_txt",
                 "exam.sty_indicator_rout_screening_obsp_yn",
                 "exam.sty_indicator_high_risk_yn",
                 "exam.sty_indicator_high_risk_brca_1_yn",  
                 "exam.sty_indicator_high_risk_brca_2_yn", 
                 "exam.sty_indicator_high_risk_brca_1_or_2_yn", 
                 "exam.sty_indicator_high_risk_at_yn", 
                 "exam.sty_indicator_high_risk_other_gene_yn", 
                 "exam.sty_indicator_high_risk_prior_high_risk_marker_yn", 
                 "exam.sty_indicator_high_risk_prior_personal_can_hist_yn", 
                 "exam.sty_indicator_high_risk_hist_of_mantle_rad_yn", 
                 "exam.sty_indicator_high_risk_fam_hist_yn", 
                 "exam.sty_indicator_add_eval_as_folup_yn", 
                 "exam.sty_indicator_folup_after_pre_exam_yn",
                 "exam.sty_indicator_pre_operative_extent_of_dis_yn",
                 "exam.sty_indicator_post_operative_margin_yn", 
                 "exam.sty_indicator_pre_neoadj_trtmnt_yn",
                 "exam.sty_indicator_prob_solv_diff_img_yn", 
                 "exam.sty_indicator_scar_vs_recurr_yn",
                 "exam.sty_indicator_folup_recommend_yn", 
                 "exam.sty_indicator_prior_2_prophy_mast_yn")
                     
    # write output query to pandas frame.
    print len(datainfo)
    print datainfo
    dinfo = pd.DataFrame(data=datainfo[0], columns=colLabels)
    # print results
    print dinfo
        
    return dinfo
        

if __name__ == '__main__':    
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__))
    print path_rootFolder
       
    # Open filename list
    print sys.argv[1]
    file_ids = open(sys.argv[1],"r")
    file_ids.seek(0)
    
    # initialize    
    casesFrame = pd.DataFrame()
    
    for k in range(1):
        line = file_ids.readline()
    print line
    
    # Create the database: the Session. 
    Session = sessionmaker()
    Session.configure(bind=engine)  # once engine is available
    session = Session() #instantiate a Session
       
    while ( line ) : 
        # Get the line:      3 column info [id  Study#  examDate]
        fileline = line.split()
        lesion_id = int(fileline[0] )
        StudyID = fileline[1] 
        examdateID = fileline[2]
       
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        print "Executing SQL connection..."
        # Format query StudyID
        if (len(StudyID) >= 4 ): fStudyID=StudyID
        if (len(StudyID) == 3 ): fStudyID='0'+StudyID
        if (len(StudyID) == 2 ): fStudyID='00'+StudyID
        if (len(StudyID) == 1 ): fStudyID='000'+StudyID
           
        # Format query reprocdateID
        redateID = datetime.date(int(examdateID[0:4]), int(examdateID[5:7]), int(examdateID[8:10]))
        dinfo = queryDatabase(session, fStudyID, redateID)
        
        ## append
        casesFrame = casesFrame.append(dinfo)
        ## continue to next case
        line = file_ids.readline()
        print "\n-------\n"+line
        
    print "\n-------\nFinalizing..."
    print casesFrame
    writer = pd.ExcelWriter('casesTesting.xlsx')
    casesFrame.to_excel(writer,'casesTesting')
    writer.save()