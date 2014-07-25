# -*- coding: utf-8 -*-
"""
Process all pipeline for a record and
Send a record to database

Created on Fri Jul 25 15:46:19 2014

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

from query_database import *
import processDicoms

from inputs_init import *
from display import *
from features_dynamic import *
from features_morphology import *
from features_texture import *
from features_T2 import *
from segment import *
import pylab      
import annot
 
from sqlalchemy import Column, Integer, String
import datetime
from mybase import myengine
import mydatabase
from sqlalchemy.orm import sessionmaker
from add_records import *


class Send(object):
    """
    USAGE:
    =============
    Send2DB = Send()
    """
    def __init__(self): 
        self.dataInfo = []
        self.queryData = Query() 
        self.load = Inputs_init()
        self.records = AddRecords()
        # Create only 1 display
        self.loadDisplay = Display()
        self.createSegment = Segment()
        
    
    def process(self, path_rootFolder, fileline, line, line_muscleVOI, lesion_id, cond, BenignNMaligNAnt, StudyID, DicomExamNumber, Lesionfile, dateID, SeriesID, finding_side, Diagnosis, T2SeriesID):
                
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        print "Executing SQL connection..."
        # Format query StudyID
        if (len(StudyID) >= 4 ): fStudyID=StudyID
        if (len(StudyID) == 3 ): fStudyID='0'+StudyID
        if (len(StudyID) == 2 ): fStudyID='00'+StudyID
        if (len(StudyID) == 1 ): fStudyID='000'+StudyID
           
        # Format query redateID
        redateID = dateID[0:4]+'-'+dateID[4:6]+'-'+dateID[6:8]
    
        # perform query        
        try:
            is_mass, colLabelsmass, is_nonmass, colLabelsnonmass = self.queryData.queryDatabase4T2(fStudyID, redateID)
            
        except Exception: 
            self.queryData.queryDatabasewNoGuiNodate(fStudyID)

        rowCase=0 
        rowCase = int(raw_input('pick row (0-n): '))
         #slice data, get only 1 record        
        dataCase = pd.Series( self.queryData.d1.loc[rowCase,:] )
        print dataCase 
        print is_mass
        print is_nonmass   
         
        ## append collection of cases
        casesFrame = pd.DataFrame(columns=self.queryData.d1.columns)
        casesFrame = casesFrame.append(dataCase) # 20    
        casesFrame['id']=fStudyID
        casesFrame.set_index('id',inplace=False) 
        
        # ask for info about lesion row data from query
        if "mass" == cond[0:-1]:
            mass_rowCase = raw_input('pick row for MASS (0-n) or x: ')
            if (mass_rowCase == 'x'):
                ## append collection of cases 
                massframe = pd.DataFrame(data=array( ['NA'] * 10))
                massframe = massframe.transpose()
                massframe.columns = ['finding.side_int','finding.size_x_double','finding.size_y_double','finding.size_z_double','finding.mri_dce_init_enh_int','finding.mri_dce_delay_enh_int','finding.curve_int','finding.mri_nonmass_dist_int', 'finding.mri_nonmass_int_enh_int','finding.t2_signal_int']
            else:
                massframe = pd.DataFrame(data=array( is_mass[int(mass_rowCase)] ))
                massframe = massframe.transpose()
                massframe.columns = list(colLabelsmass)
            
        if "nonmass" == cond[0:-1]:
            nonmass_rowCase = raw_input('pick row for NONMASS (0-n) or x: ')
            if (nonmass_rowCase == 'x'):
                ## append collection of cases 
                nonmassframe = pd.DataFrame(data=array( ['NA'] * 10))
                nonmassframe = nonmassframe.transpose()
                nonmassframe.columns = ['finding.side_int','finding.size_x_double','finding.size_y_double','finding.size_z_double','finding.mri_dce_init_enh_int','finding.mri_dce_delay_enh_int','finding.curve_int','finding.mri_nonmass_dist_int', 'finding.mri_nonmass_int_enh_int','finding.t2_signal_int']
            else:
                nonmassframe = pd.DataFrame(data=array( is_nonmass[int(nonmass_rowCase)] ))
                nonmassframe = nonmassframe.transpose()
                nonmassframe.columns = list(colLabelsnonmass)
        
        #slice data, get only 1 record    
        if "mass" == cond[0:-1]:
            dataCase = massframe 
            img_folder = 'Z:/Cristina/MassNonmass/mass/'
            casesFrame['finding.mri_mass_yn'] = True
            casesFrame['finding.mri_nonmass_yn'] = False
            
        if "nonmass" == cond[0:-1]: 
            dataCase = nonmassframe 
            img_folder = 'Z:/Cristina/MassNonmass/nonmass/'
            casesFrame['finding.mri_nonmass_yn'] = True
            casesFrame['finding.mri_mass_yn'] = False
        
        # analyze T2 series            
        if (T2SeriesID == 'NONE'):
             # ask for which series to load
            [abspath_ExamID, eID, SeriesIDall, studyFolder, dicomInfo] = processDicoms.get_series(StudyID, img_folder)
            
            print "\n----------------------------------------------------------"
            choseSerie = raw_input('Enter n T2 Series to load (0-n), or x if NO T2w sequence and pass: ')
            if (choseSerie == 'x'):
                ## append collection of cases 
                casesFrame['T2SeriesID']='NONE'
                casesFrame['cond']=cond
                casesFrame['id']=fStudyID
                casesFrame['DicomExamNumber']=DicomExamNumber
                casesFrame['Lesionfile']=Lesionfile
                casesFrame['BenignNMaligNAnt']=BenignNMaligNAnt
                casesFrame['Diagnosis']=Diagnosis
            else:    
                T2SeriesID = SeriesIDall[int(choseSerie)]
                
                ## append collection of cases 
                casesFrame['T2SeriesID']=T2SeriesID
                casesFrame['cond']=cond
                casesFrame['id']=fStudyID
                casesFrame['DicomExamNumber']=DicomExamNumber
                casesFrame['Lesionfile']=Lesionfile
                casesFrame['BenignNMaligNAnt']=BenignNMaligNAnt
                casesFrame['Diagnosis']=Diagnosis
        
                path_T2Series = img_folder+StudyID+'/'+eID+'/'+T2SeriesID
                print "\nPath to T2 location: %s" %  path_T2Series
        else:
            path_T2Series = img_folder+StudyID+'/'+DicomExamNumber+'/'+T2SeriesID
            print "\nPath to T2 location: %s" %  path_T2Series
        
        #############################                  
        ###### 2) Check segmentation accuracy with annotations
        #############################
        ###### Start by Loading 
        print "Start by loading volumes..."
        data_loc='Z:\Cristina\MassNonmass'+os.sep+cond[:-1]
        
        [series_path, phases_series, lesionID_path] = self.load.readVolumes(data_loc, StudyID, DicomExamNumber, SeriesID, Lesionfile)
        print "Path to series location: %s" % series_path 
        print "List of pre and post contrast volume names: %s" % phases_series
        print "Path to lesion segmentation: %s" % lesionID_path
        
        print "\n Load Segmentation..."
        lesion3D = self.load.loadSegmentation(lesionID_path)
        print "Data Structure: %s" % lesion3D.GetClassName()
        print "Number of points: %d" % int(lesion3D.GetNumberOfPoints())
        print "Number of cells: %d" % int(lesion3D.GetNumberOfCells())
        
        #############################
        # 4) Parse annotations (display and pick corresponding to lesion)
        #############################
        print "\n Visualize volumes..."
           
        # extract annotation if any    
        if len(fileline) > 13:
            if fileline[13] !='NA' and fileline[13] !='[]'  :
                annots_all = line[line.find('[')+1:line.find(']')] 
                annots = True
                list_annots=[]
                # iterate throuhg
                while(annots):
                    one_annot = annots_all[annots_all.find('{')+1:annots_all.find('}')]
                    list_annots.append(one_annot)
                    annots_all = annots_all[annots_all.find('}')+1:]
                    if '{' not in annots_all:   annots = False
                
                print "\nLoading annotations..." 
                annots_dict_list = self.loadDisplay.extract_annot(list_annots)
                print "\nDisplay annotations:" 
                self.loadDisplay.display_annot(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, annots_dict_list)
                    
        else:
            list_annots = 'NA'
            print "\n####################"
            print "No Annotations"
            print "####################"
                             
        self.loadDisplay.addSegment(lesion3D, (0,1,0))
        self.createSegment.saveSegmentation(path_rootFolder+os.sep+'segmentations', lesion3D, lesionfilename=StudyID+'_'+DicomExamNumber+'_'+Lesionfile+'.vtk') 
        self.loadDisplay.visualize(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, sub=True, postS=4, interact=True)
       
        chooseAnnot = int(raw_input('\n Enter # corresponding to Lesion Annotation or 0 to skip: ') )
        if chooseAnnot != 0:
            casesFrame['LesionAnnot']= str(annots_dict_list[chooseAnnot-1])
            annot_attrib = annots_dict_list[chooseAnnot-1]
            pi = annot_attrib['pi_2display']
            pf = annot_attrib['pf_2display'] 
            
            #############################
            ###### Compare manual marker distance with auto segmentation length for validation
            #############################
            eu_dist_mkers = float( sqrt( (pi[0]-pf[0])**2 + (pi[1]-pf[1])**2 + (pi[2]-pf[2])**2 ) )           
            print "eu_dist_mkers: " 
            print eu_dist_mkers
            
            axis_lengths = self.loadDisplay.extract_segment_dims(lesion3D)
            eu_dist_seg =  float(sqrt( axis_lengths[0] + axis_lengths[1]))  # only measure x-y euclidian distance betweeen extreme points
            print "eu_dist_seg : " 
            print eu_dist_seg 
       
        else:
            annot_attrib=[]
            casesFrame['LesionAnnot']= "[]"
                
        #############################
        # 4) Manually modify Segmentation of lesion. Comment out if not needed ( define seededlesion3D = lesion3D  )
        #############################
        chgSeg = int(raw_input('\n Enter 1 to modify segmentation, 0 to skip: ') )
        if chgSeg == 1:
            #  Get z slice
            LesionZslice = self.loadDisplay.zImagePlaneWidget.GetSliceIndex()
            
            print "\n Displaying picker for lesion segmentation"
            seeds = self.loadDisplay.display_pick(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, 4, self.LesionZslice)
            
            seededlesion3D = self.createSegment.segmentFromSeeds(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, seeds, self.loadDisplay.iren1, self.loadDisplay.xImagePlaneWidget, self.loadDisplay.yImagePlaneWidget,  self.loadDisplay.zImagePlaneWidget)
            self.loadDisplay.addSegment(seededlesion3D, (0,0,1))
            self.loadDisplay.picker.RemoveAllObservers()
            
            # save it to file	             
            createSegment.saveSegmentation(lesionID_path, seededlesion3D) 
            createSegment.saveSegmentation(path_rootFolder+os.sep+'segmentations', seededlesion3D, lesionfilename=StudyID+'_'+DicomExamNumber+'_'+Lesionfile+'.vtk') 
            lesion3D = seededlesion3D
            
            axis_lengths = loadDisplay.extract_segment_dims(lesion3D)
            print axis_lengths
            eu_dist_seg = float( sqrt( axis_lengths[0] + axis_lengths[1])) # only measure x-y euclidian distance betweeen extreme points
            print "eu_dist_seg : " 
            print eu_dist_seg 
            
            loadDisplay.renderer1.RemoveActor(loadDisplay.actor_mesh)
                      
        #############################
        ###### Extract Dynamic features
        #############################
        print "\n Extract Dynamic contour features..."
        loadDynamic = Dynamic()
        dyn_contour = loadDynamic.extractfeatures_contour(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
        print "\n=========================================="
        print dyn_contour
                
        print "\n Extract Dynamic inside features..."
        dyn_inside = loadDynamic.extractfeatures_inside(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
        print dyn_inside
        print "\n=========================================="
 
        pylab.close('all') 
         
        #############################
        ###### Extract Morphology features
        #############################
        print "\n Extract Morphology features..."
        loadMorphology = Morphology()
        morphofeatures = loadMorphology.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D)
        print "\n=========================================="
        print morphofeatures
        print "\n=========================================="

        pylab.close('all') 
                
        #############################        
        ###### Extract Texture features
        #############################
        print "\n Extract Texture features..."
        loadTexture = Texture()
        texturefeatures = loadTexture.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, lesion3D, loadMorphology.VOI_efect_diameter, loadMorphology.lesion_centroid )
        print "\n=========================================="
        print texturefeatures
        print "\n=========================================="

        pylab.close('all')  
                         
        #############################
        # Extract Lesion and Muscle Major pectoralies signal                                   
        ############################# 
        if T2SeriesID != 'NONE':     
            #############################        
            ###### Extract T2 features, Process T2 and visualize
            #############################
            ###### Start by Loading 
            print "Start by loading T2 volume..."       
            load.readT2(path_T2Series)
            
            print "\n Visualize addT2visualize ..."
            loadDisplay.addT2visualize(load.T2Images, load.T2image_pos_pat, load.T2image_ori_pat, load.T2dims, load.T2spacing, interact=True)

            T2 = features_T2()  
            
            if line_muscleVOI.split()[0] != "NA" and line_muscleVOI.split()[0] != '[]':
                l = line_muscleVOI[line_muscleVOI.find('[')+1:-2].split(",")
                bounds_muscleSI = [float(l[0]), float(l[1]), float(l[2]), float(l[3]), float(l[4]), float(l[5]) ]
                print "\n bounds_muscleSI from file:"
                print bounds_muscleSI
                       
                # instead of extract_muscleSI use load_muscleSI from file
                [T2_muscleSI, muscle_scalar_range]  = T2.load_muscleSI(load.T2Images, load.T2image_pos_pat, load.T2image_ori_pat, loadDisplay.origin, loadDisplay.iren1, loadDisplay.renderer1, loadDisplay.picker, loadDisplay.xImagePlaneWidget, loadDisplay.yImagePlaneWidget, loadDisplay.zImagePlaneWidget, bounds_muscleSI)
                print "ave. T2_muscleSI: %d" % mean(T2_muscleSI)
    
            else:
                # Do extract_muscleSI 
                [T2_muscleSI, muscle_scalar_range]  = T2.extract_muscleSI(load.T2Images, load.T2image_pos_pat, load.T2image_ori_pat, loadDisplay.origin, loadDisplay.iren1, loadDisplay.renderer1, loadDisplay.picker, loadDisplay.xImagePlaneWidget, loadDisplay.yImagePlaneWidget, loadDisplay.zImagePlaneWidget)
                print "ave. T2_muscleSI: %d" % mean(T2_muscleSI)
                
            loadDisplay.iren1.Start()
            
            # Do extract_lesionSI       
            [T2_lesionSI, lesion_scalar_range]  = T2.extract_lesionSI(load.T2Images, lesion3D, load.T2image_pos_pat, load.T2image_ori_pat)
            print "ave. T2_lesionSI: %d" % mean(T2_lesionSI)
            
            LMSIR = mean(T2_lesionSI)/mean(T2_muscleSI)
            print "LMSIR: %d" % LMSIR
                    
            #############################
            # Extract morphological and margin features from T2                                   
            #############################
            print "\n Extract T2 Morphology features..."
            morphoT2features = T2.extractT2morphology(load.T2Images, lesion3D, load.T2image_pos_pat, load.T2image_ori_pat)
            print "\n=========================================="
            print morphoT2features
            print "\n Extract T2 Texture features..."
            textureT2features = T2.extractT2texture(load.T2Images, lesion3D, load.T2image_pos_pat, load.T2image_ori_pat)
            print textureT2features
            print "\n=========================================="
            

            
        # finally transform centroid world coords to ijk indexes
        im_pt = [0,0,0]
        ijk = [0,0,0]
        pco = [0,0,0]
        pixId_sliceloc = loadDisplay.transformed_image.FindPoint(loadDisplay.lesion_centroid)
        loadDisplay.transformed_image.GetPoint(pixId_sliceloc, im_pt) 
        io = loadDisplay.transformed_image.ComputeStructuredCoordinates( im_pt, ijk, pco)
        if io:
            lesion_centroid_ijk = ijk
            
        # deal with closing windows, plots, renders, actors
        pylab.close('all')
        loadDisplay.renderer1.FastDelete()
        loadDisplay.renWin1.Finalize()
        loadDisplay.iren1.TerminateApp()
                
        #############################
        ###### Finish tidying up a"nd save to file
        ## append collection of cases
        #############################  
        print "\n Adding record case to DB..."
        self.records.lesion_2DB(Lesionfile, fStudyID, DicomExamNumber, str(casesFrame['exam.a_number_txt'][rowCase]), datetime.date(int(dateID[0:4]), int(dateID[4:6]), int(dateID[6:8])), str(casesFrame['exam.mri_cad_status_txt'][rowCase]), 
                           str(casesFrame['cad.latest_mutation'][rowCase]), casesFrame['finding.mri_mass_yn'][rowCase], casesFrame['finding.mri_nonmass_yn'][rowCase], finding_side, str(casesFrame['proc.pt_procedure_id'][rowCase]), 
                            casesFrame['proc.proc_dt_datetime'][rowCase], str(casesFrame['proc.proc_side_int'][rowCase]), str(casesFrame['proc.proc_source_int'][rowCase]),  str(casesFrame['proc.proc_guid_int'][rowCase]), 
                            str(casesFrame['proc.proc_tp_int'][rowCase]), str(casesFrame['exam.comment_txt'][rowCase]), str(casesFrame['proc.original_report_txt'][rowCase]), str(dataCase['finding.curve_int'][0]), str(dataCase['finding.mri_dce_init_enh_int'][0]), str(dataCase['finding.mri_dce_delay_enh_int'][0]), cond,  Diagnosis)
        
        if "mass" == cond[0:-1]:
            self.records.mass_2DB(lesion_id, str(BenignNMaligNAnt), SeriesID, T2SeriesID, dataCase['finding.mammo_n_mri_mass_shape_int'][0], dataCase['finding.mri_mass_margin_int'][0] )
            rowid = mass_rowCase

        if "nonmass" == cond[0:-1]: 
            self.records.nonmass_2DB(lesion_id, str(BenignNMaligNAnt), SeriesID, T2SeriesID, dataCase['finding.mri_nonmass_dist_int'][0], dataCase['finding.mri_nonmass_int_enh_int'][0])
            rowid = nonmass_rowCase
            
        # send features
        # Dynamic
        self.records.dyn_records_2DB(lesion_id, dyn_inside['A.inside'], dyn_inside['alpha.inside'], dyn_inside['beta.inside'], dyn_inside['iAUC1.inside'], dyn_inside['Slope_ini.inside'], dyn_inside['Tpeak.inside'], dyn_inside['Kpeak.inside'], dyn_inside['SER.inside'], dyn_inside['maxCr.inside'], dyn_inside['peakCr.inside'], dyn_inside['UptakeRate.inside'], dyn_inside['washoutRate.inside'], dyn_inside['maxVr.inside'], dyn_inside['peakVr.inside'], dyn_inside['Vr_increasingRate.inside'], dyn_inside['Vr_decreasingRate.inside'], dyn_inside['Vr_post_1.inside'],
                               dyn_contour['A.contour'], dyn_contour['alpha.contour'], dyn_contour['beta.contour'], dyn_contour['iAUC1.contour'], dyn_contour['Slope_ini.contour'], dyn_contour['Tpeak.contour'], dyn_contour['Kpeak.contour'], dyn_contour['SER.contour'], dyn_contour['maxCr.contour'], dyn_contour['peakCr.contour'], dyn_contour['UptakeRate.contour'], dyn_contour['washoutRate.contour'], dyn_contour['maxVr.contour'], dyn_contour['peakVr.contour'], dyn_contour['Vr_increasingRate.contour'], dyn_contour['Vr_decreasingRate.contour'], dyn_contour['Vr_post_1.contour'] )
        
        # Morphology
        self.records.morpho_records_2DB(lesion_id, morphofeatures['min_F_r_i'], morphofeatures['max_F_r_i'], morphofeatures['mean_F_r_i'], morphofeatures['var_F_r_i'], morphofeatures['skew_F_r_i'], morphofeatures['kurt_F_r_i'], morphofeatures['iMax_Variance_uptake'], 
                                                  morphofeatures['iiMin_change_Variance_uptake'], morphofeatures['iiiMax_Margin_Gradient'], morphofeatures['k_Max_Margin_Grad'], morphofeatures['ivVariance'], morphofeatures['circularity'], morphofeatures['irregularity'], morphofeatures['edge_sharp_mean'],
                                                  morphofeatures['edge_sharp_std'], morphofeatures['max_RGH_mean'], morphofeatures['max_RGH_mean_k'], morphofeatures['max_RGH_var'], morphofeatures['max_RGH_var_k'] )
        # Texture
        self.records.texture_records_2DB(lesion_id, texturefeatures['texture_contrast_zero'], texturefeatures['texture_contrast_quarterRad'], texturefeatures['texture_contrast_halfRad'], texturefeatures['texture_contrast_threeQuaRad'], 
                                                  texturefeatures['texture_homogeneity_zero'], texturefeatures['texture_homogeneity_quarterRad'], texturefeatures['texture_homogeneity_halfRad'], texturefeatures['texture_homogeneity_threeQuaRad'], 
                                                  texturefeatures['texture_dissimilarity_zero'], texturefeatures['texture_dissimilarity_quarterRad'], texturefeatures['texture_dissimilarity_halfRad'], texturefeatures['texture_dissimilarity_threeQuaRad'], 
                                                  texturefeatures['texture_correlation_zero'], texturefeatures['texture_correlation_quarterRad'], texturefeatures['texture_correlation_halfRad'], texturefeatures['texture_correlation_threeQuaRad'], 
                                                  texturefeatures['texture_ASM_zero'], texturefeatures['texture_ASM_quarterRad'], texturefeatures['texture_ASM_halfRad'], texturefeatures['texture_ASM_threeQuaRad'], 
                                                  texturefeatures['texture_energy_zero'], texturefeatures['texture_energy_quarterRad'], texturefeatures['texture_energy_halfRad'], texturefeatures['texture_energy_threeQuaRad'] )
        # T2 relative signal, morphology and texture
        if T2SeriesID != 'NONE':                                                       
            self.records.t2_records_2DB(lesion_id, dataCase['finding.t2_signal_int'][0], str(list(load.T2dims)), str(list(load.T2spacing)), str(load.T2fatsat), mean(T2_muscleSI), std(T2_muscleSI), str(muscle_scalar_range), str(T2.bounds_muscleSI), mean(T2_lesionSI), std(T2_lesionSI), str(lesion_scalar_range), LMSIR, 
                                            morphoT2features['T2min_F_r_i'], morphoT2features['T2max_F_r_i'], morphoT2features['T2mean_F_r_i'], morphoT2features['T2var_F_r_i'], morphoT2features['T2skew_F_r_i'], morphoT2features['T2kurt_F_r_i'], morphoT2features['T2grad_margin'], morphoT2features['T2grad_margin_var'], morphoT2features['T2RGH_mean'], morphoT2features['T2RGH_var'], 
                                            textureT2features['T2texture_contrast_zero'], textureT2features['T2texture_contrast_quarterRad'], textureT2features['T2texture_contrast_halfRad'], textureT2features['T2texture_contrast_threeQuaRad'], 
                                            textureT2features['T2texture_homogeneity_zero'], textureT2features['T2texture_homogeneity_quarterRad'], textureT2features['T2texture_homogeneity_halfRad'], textureT2features['T2texture_homogeneity_threeQuaRad'], 
                                            textureT2features['T2texture_dissimilarity_zero'], textureT2features['T2texture_dissimilarity_quarterRad'], textureT2features['T2texture_dissimilarity_halfRad'], textureT2features['T2texture_dissimilarity_threeQuaRad'], 
                                            textureT2features['T2texture_correlation_zero'], textureT2features['T2texture_correlation_quarterRad'], textureT2features['T2texture_correlation_halfRad'], textureT2features['T2texture_correlation_threeQuaRad'], 
                                            textureT2features['T2texture_ASM_zero'], textureT2features['T2texture_ASM_quarterRad'], textureT2features['T2texture_ASM_halfRad'], textureT2features['T2texture_ASM_threeQuaRad'], 
                                            textureT2features['T2texture_energy_zero'], textureT2features['T2texture_energy_quarterRad'], textureT2features['T2texture_energy_halfRad'], textureT2features['T2texture_energy_threeQuaRad'])
        # Send annotation if any
        if annot_attrib:
            self.records.annot_records_2DB(lesion_id, annot_attrib['AccessionNumber'], annot_attrib['SeriesDate'], annot_attrib['SeriesNumber'], annot_attrib['SliceLocation'], annot_attrib['SeriesDescription'], annot_attrib['PatientID'], annot_attrib['StudyID'], annot_attrib['SeriesInstanceUID'], annot_attrib['note'], annot_attrib['xi'], annot_attrib['yi'], annot_attrib['xf'], annot_attrib['yf'], 
                                                    str(annot_attrib['pi_ijk']), str(annot_attrib['pi_2display']), str(annot_attrib['pf_ijk']), str(annot_attrib['pf_2display']),
                                                    eu_dist_mkers, eu_dist_seg)
        
        # SEgmentation details
        self.records.segment_records_2DB(lesion_id, loadDisplay.lesion_bounds[0], loadDisplay.lesion_bounds[1], loadDisplay.lesion_bounds[2], loadDisplay.lesion_bounds[3], loadDisplay.lesion_bounds[4], loadDisplay.lesion_bounds[5],
                                                    loadDisplay.no_pts_segm, loadDisplay.VOI_vol, loadDisplay.VOI_surface, loadDisplay.VOI_efect_diameter, str(list(loadDisplay.lesion_centroid)), str(lesion_centroid_ijk))
                                                    
        return