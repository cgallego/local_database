# -*- coding: utf-8 -*-
"""
USAGE: 
=============
from add_records import *
record = AddRecords()
record.lesion_2DB(image, image_pos_pat, image_ori_pat)
record.addSegment(lesion3D)
record.subImage(Images2Sub, timep)                  
record.visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Methods:
=============
dicomTransform(image, image_pos_pat, image_ori_pat)
addSegment(lesion3D)
subImage(Images2Sub, timep)                  
visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Instance Attributes:
===============

Created on Wed Jul 23 09:32:26 2014

@ author (C) Cristina Gallego, University of Toronto
--------------------------------------------------------------------
 """
import os, os.path
import sys
import string
from sys import argv, stderr, exit
import vtk
from numpy import *

from sqlalchemy.orm import sessionmaker
from mybase import myengine
import mydatabase

class AddRecords(object):
    """
    USAGE:
    =============
    record = AddRecords()
    """
    def __init__(self): 
        """ initialize database session """           
        #  create a top level Session configuration which can then be used throughout
        # Create the Session
        self.Session = sessionmaker()
        self.Session.configure(bind=myengine)  # once engine is available
        
    def __call__(self):       
        """ Turn Class into a callable object """
        AddRecords() 

    def lesion_2DB(self, lesionfile, cad_id, dicom_no, accession_no, exam_date, cad_status, mutation, mass_yn, nonmass_yn, finding_side, proc_id, proc_date, proc_side, proc_source, proc_guid, proc_type, lesion_comments, original_report, curve_int, dce_init, dce_delay, label, diagnosis):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        lesion_info = mydatabase.Lesion_record(lesionfile, cad_id, dicom_no, accession_no, exam_date, cad_status, mutation, mass_yn, nonmass_yn, finding_side, proc_id, proc_date, proc_side, proc_source, proc_guid, proc_type, lesion_comments, original_report, curve_int, dce_init, dce_delay, label, diagnosis)
        self.session.add(lesion_info)
        print self.session.query(mydatabase.Lesion_record).first()
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return
        
    def mass_2DB(self, lesion_id, BenignNMaligNAnt, DynSeries_id, T2Series_id, mri_mass_shape, mri_mass_margin):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        mass_record = mydatabase.Mass_record(lesion_id, BenignNMaligNAnt, DynSeries_id, T2Series_id, mri_mass_shape, mri_mass_margin)
        self.session.add(mass_record)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return        
            
    def nonmass_2DB(self, lesion_id, BenignNMaligNAnt, DynSeries_id, T2Series_id, mri_nonmass_dist, mri_nonmass_int_enh):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        nonmass_record = mydatabase.Nonmass_record(lesion_id, BenignNMaligNAnt, DynSeries_id, T2Series_id, mri_nonmass_dist, mri_nonmass_int_enh)
        self.session.add(nonmass_record)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return
        
    def dyn_records_2DB(self, lesion_id, A_inside, alpha_inside, beta_inside, iAUC1_inside, Slope_ini_inside, Tpeak_inside,
                 Kpeak_inside, SER_inside, maxCr_inside, peakCr_inside, UptakeRate_inside, washoutRate_inside, maxVr_inside,
                 peakVr_inside, Vr_increasingRate_inside, Vr_decreasingRate_inside, Vr_post_1_inside,
                 A_countor, alpha_countor, beta_countor, iAUC1_countor, Slope_ini_countor, Tpeak_countor,
                 Kpeak_countor, SER_countor, maxCr_countor, peakCr_countor, UptakeRate_countor, washoutRate_countor, maxVr_countor,
                 peakVr_countor, Vr_increasingRate_countor, Vr_decreasingRate_countor, Vr_post_1_countor):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        dyn_records = mydatabase.Dynamic_features(lesion_id, A_inside, alpha_inside, beta_inside, iAUC1_inside, Slope_ini_inside, Tpeak_inside,
                 Kpeak_inside, SER_inside, maxCr_inside, peakCr_inside, UptakeRate_inside, washoutRate_inside, maxVr_inside,
                 peakVr_inside, Vr_increasingRate_inside, Vr_decreasingRate_inside, Vr_post_1_inside,
                 A_countor, alpha_countor, beta_countor, iAUC1_countor, Slope_ini_countor, Tpeak_countor,
                 Kpeak_countor, SER_countor, maxCr_countor, peakCr_countor, UptakeRate_countor, washoutRate_countor, maxVr_countor,
                 peakVr_countor, Vr_increasingRate_countor, Vr_decreasingRate_countor, Vr_post_1_countor)
        self.session.add(dyn_records)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return

    def morpho_records_2DB(self, lesion_id, min_F_r_i, max_F_r_i, mean_F_r_i, var_F_r_i, skew_F_r_i, kurt_F_r_i,
                 iMax_Variance_uptake, iiMin_change_Variance_uptake, iiiMax_Margin_Gradient, k_Max_Margin_Grad,
                 ivVariance, circularity, irregularity, edge_sharp_mean, edge_sharp_std, max_RGH_mean, max_RGH_mean_k, max_RGH_var, max_RGH_var_k):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        morp_records = mydatabase.Morpho_features(lesion_id, min_F_r_i, max_F_r_i, mean_F_r_i, var_F_r_i, skew_F_r_i, kurt_F_r_i,
                 iMax_Variance_uptake, iiMin_change_Variance_uptake, iiiMax_Margin_Gradient, k_Max_Margin_Grad,
                 ivVariance, circularity, irregularity, edge_sharp_mean, edge_sharp_std, max_RGH_mean, max_RGH_mean_k, max_RGH_var, max_RGH_var_k)
        self.session.add(morp_records)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return
        
    
    def texture_records_2DB(self, lesion_id, texture_contrast_zero, texture_contrast_quarterRad, texture_contrast_halfRad, texture_contrast_threeQuaRad,  
                 texture_homogeneity_zero, texture_homogeneity_quarterRad, texture_homogeneity_halfRad, texture_homogeneity_threeQuaRad,
                 texture_dissimilarity_zero, texture_dissimilarity_quarterRad, texture_dissimilarity_halfRad, texture_dissimilarity_threeQuaRad,
                 texture_correlation_zero, texture_correlation_quarterRad, texture_correlation_halfRad, texture_correlation_threeQuaRad,
                 texture_ASM_zero, texture_ASM_quarterRad, texture_ASM_halfRad, texture_ASM_threeQuaRad,
                 texture_energy_zero, texture_energy_quarterRad, texture_energy_halfRad, texture_energy_threeQuaRad):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        tex_records = mydatabase.Texture_features(lesion_id, texture_contrast_zero, texture_contrast_quarterRad, texture_contrast_halfRad, texture_contrast_threeQuaRad,  
                 texture_homogeneity_zero, texture_homogeneity_quarterRad, texture_homogeneity_halfRad, texture_homogeneity_threeQuaRad,
                 texture_dissimilarity_zero, texture_dissimilarity_quarterRad, texture_dissimilarity_halfRad, texture_dissimilarity_threeQuaRad,
                 texture_correlation_zero, texture_correlation_quarterRad, texture_correlation_halfRad, texture_correlation_threeQuaRad,
                 texture_ASM_zero, texture_ASM_quarterRad, texture_ASM_halfRad, texture_ASM_threeQuaRad,
                 texture_energy_zero, texture_energy_quarterRad, texture_energy_halfRad, texture_energy_threeQuaRad)
        self.session.add(tex_records)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return
        

    def t2_records_2DB(self, lesion_id, find_t2_signal_int, T2dims, T2spacing, T2fatsat,  
                 T2_muscleSI, T2_muscleSIstd, muscle_scalar_range, bounds_muscleSI,
                 T2_lesionSI, T2_lesionSIstd, lesion_scalar_range, LMSIR,
                 T2min_F_r_i, T2max_F_r_i, T2mean_F_r_i, T2var_F_r_i, T2skew_F_r_i, T2kurt_F_r_i,
                 T2grad_margin, T2grad_margin_var, T2RGH_mean, T2RGH_var,
                 T2texture_contrast_zero, T2texture_contrast_quarterRad, T2texture_contrast_halfRad, T2texture_contrast_threeQuaRad,  
                 T2texture_homogeneity_zero, T2texture_homogeneity_quarterRad, T2texture_homogeneity_halfRad, T2texture_homogeneity_threeQuaRad,
                 T2texture_dissimilarity_zero, T2texture_dissimilarity_quarterRad, T2texture_dissimilarity_halfRad, T2texture_dissimilarity_threeQuaRad,
                 T2texture_correlation_zero, T2texture_correlation_quarterRad, T2texture_correlation_halfRad, T2texture_correlation_threeQuaRad,
                 T2texture_ASM_zero, T2texture_ASM_quarterRad, T2texture_ASM_halfRad, T2texture_ASM_threeQuaRad,
                 T2texture_energy_zero, T2texture_energy_quarterRad, T2texture_energy_halfRad, T2texture_energy_threeQuaRad):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        t2_records = mydatabase.T2_features(lesion_id, find_t2_signal_int, T2dims, T2spacing, T2fatsat,  
                 T2_muscleSI, T2_muscleSIstd, muscle_scalar_range, bounds_muscleSI,
                 T2_lesionSI, T2_lesionSIstd, lesion_scalar_range, LMSIR,
                 T2min_F_r_i, T2max_F_r_i, T2mean_F_r_i, T2var_F_r_i, T2skew_F_r_i, T2kurt_F_r_i,
                 T2grad_margin, T2grad_margin_var, T2RGH_mean, T2RGH_var,
                 T2texture_contrast_zero, T2texture_contrast_quarterRad, T2texture_contrast_halfRad, T2texture_contrast_threeQuaRad,  
                 T2texture_homogeneity_zero, T2texture_homogeneity_quarterRad, T2texture_homogeneity_halfRad, T2texture_homogeneity_threeQuaRad,
                 T2texture_dissimilarity_zero, T2texture_dissimilarity_quarterRad, T2texture_dissimilarity_halfRad, T2texture_dissimilarity_threeQuaRad,
                 T2texture_correlation_zero, T2texture_correlation_quarterRad, T2texture_correlation_halfRad, T2texture_correlation_threeQuaRad,
                 T2texture_ASM_zero, T2texture_ASM_quarterRad, T2texture_ASM_halfRad, T2texture_ASM_threeQuaRad,
                 T2texture_energy_zero, T2texture_energy_quarterRad, T2texture_energy_halfRad, T2texture_energy_threeQuaRad)
        self.session.add(t2_records)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return
        
    def annot_records_2DB(self,  lesion_id, AccessionNumber, SeriesDate, SeriesNumber, SliceLocation, SeriesDescription,
                 PatientID, StudyID, SeriesInstanceUID, note, xi_coord, yi_coord, xf_coord, yf_coord,
                 pi_ijk, pi_2display, pf_ijk, pf_2display, eu_dist_mkers, eu_dist_seg):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        annot_records = mydatabase.Annot_record( lesion_id, AccessionNumber, SeriesDate, SeriesNumber, SliceLocation, SeriesDescription,
                 PatientID, StudyID, SeriesInstanceUID, note, xi_coord, yi_coord, xf_coord, yf_coord,
                 pi_ijk, pi_2display, pf_ijk, pf_2display, eu_dist_mkers, eu_dist_seg)
        self.session.add(annot_records)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return
        

    def segment_records_2DB(self, lesion_id, segm_xmin, segm_xmax, segm_ymin, segm_ymax, segm_zmin, segm_zmax,
                 no_pts, voi_vol, voi_surface, VOI_efect_diameter, lesion_centroid_world, lesion_centroid_ijk):        
        
        self.session = self.Session() #instantiate a Session
        # Send to database lesion info
        segment_records = mydatabase.Segment_record(lesion_id, segm_xmin, segm_xmax, segm_ymin, segm_ymax, segm_zmin, segm_zmax,
                 no_pts, voi_vol, voi_surface, VOI_efect_diameter, lesion_centroid_world, lesion_centroid_ijk)
        self.session.add(segment_records)
        
        # Finally send records to database
        try:
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()
            
        return     