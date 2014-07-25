# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:48:07 2014

@ author (C) Cristina Gallego, University of Toronto
"""
import sys, os
import string
import datetime
from numpy import *


import datetime
from mybase import Base, myengine
import mydatabase

from sqlalchemy.orm import sessionmaker

import wx
import wxTableBase

#!/usr/bin/env python
class Query(object):
    """
    USAGE:
    =============
    localdata = QuerymyDatabase()
    """
    def __call__(self):       
        """ Turn Class into a callable object """
        QuerymyDatabase()
        
        
    def __init__(self): 
        """ initialize QueryDatabase """
        
        
        
                              
    def QuerymyDatabase(self):
        """
        run : Query by StudyID/AccesionN pair study to local folder
        
        Inputs
        ======
        StudyID : (int)    CAD StudyID
        redateID : (int)  CAD StudyID Data of exam (format yyyy-mm-dd)
        
        Output
        ======
        """               
        # Create the database: the Session. 
        self.Session = sessionmaker()
        self.Session.configure(bind=myengine)  # once engine is available
        session = self.Session() #instantiate a Session
        
        # Create first display
        """ Creates Table grid and Query output display Cad_Container"""
        self.app = wx.App(False)
        self.display = wxTableBase.Container(self.app)
        
       # for cad_case in session.query(Cad_record).order_by(Cad_record.pt_id): 
       #     print cad_case.pt_id, cad_case.cad_pt_no_txt, cad_case.latest_mutation_status_int    
        
        datainfo = []
        for lesion, exam, finding, proc, patho in session.query(database.Cad_record, database.Exam_record, database.Exam_Finding, database.Procedure, database.Pathology).\
                     filter(database.Cad_record.pt_id==database.Exam_record.pt_id).\
                     filter(database.Exam_record.pt_exam_id==database.Exam_Finding.pt_exam_id).\
                     filter(database.Exam_record.pt_id==database.Procedure.pt_id).\
                     filter(database.Procedure.pt_procedure_id==database.Pathology.pt_procedure_id).\
                     filter(database.Cad_record.cad_pt_no_txt == str(StudyID)).\
                     filter(database.Exam_record.exam_dt_datetime == str(redateID)).all():
                         
           # print results
           if not cad:
               print "cad is empty"
           if not exam:
               print "exam is empty"
           if not finding:
               print "finding is empty"
           if not proc:
               print "proc is empty"
           if not patho:
               print "patho is empty"
                   
           datainfo.append([cad.cad_pt_no_txt, cad.latest_mutation_status_int,
              exam.exam_dt_datetime, exam.a_number_txt, exam.mri_cad_status_txt, exam.comment_txt,
              finding.mri_mass_yn, finding.mri_nonmass_yn, finding.mri_foci_yn,
              proc.pt_procedure_id, proc.proc_dt_datetime, proc.proc_side_int, proc.proc_source_int, proc.proc_guid_int, proc.proc_tp_int, proc.original_report_txt])
           
           
           

