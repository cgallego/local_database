# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:09:01 2014

@author: Cristina Gallego
"""
import sys, os
import string
import datetime
import numpy as np

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref

from mybase import Base

#  created a Cad_record mapping 
class Lesion_record(Base):
    """Base for Exam_record class using Declarative. for table lesion
    attributes:
        self.lesion_id = lesion_id
        self.cad_pt_no_txt = cad_id
        self.exam_img_dicom_txt = dicom_no
        self.exam_a_number_txt = accession_no
        self.exam_dt_datetime = exam_date
        self.exam_mri_cad_status_txt = cad_status
        self.cad_latest_mutation_status_int = mutation
        self.exam_find_mri_mass_yn = mass_yn
        self.exam_find_mri_nonmass_yn = nonmass_yn
        self.exam_find_side_int = finding_side
        self.proc_pt_procedure_id = proc_id
        self.proc_proc_dt_datetime = proc_date
        self.proc_proc_side_int = proc_side
        self.proc_proc_source_int = proc_source
        self.proc_proc_guid_int = proc_guid
        self.proc_proc_tp_int = proc_type
        self.proc_lesion_comments_txt = lesion_comments
        self.proc_original_report_txt = original_report
        self.find_curve_int = curve_int
        self.find_mri_dce_init_enh_int = dce_init
        self.find_mri_dce_delay_enh_int = dce_delay
        self.lesion_label = label
        self.lesion_diagnosis = diagnosis    
    """
    __tablename__ = 'lesion'
    __table_args__ = {'autoload':True}
    lesion_id = Column(Integer, primary_key=True)
    mass_lesion = relationship("Mass_record", backref=backref('lesion', order_by=lesion_id))
    nonmass_lesion = relationship("Nonmass_record", backref=backref('lesion', order_by=lesion_id))
    f_dynamic = relationship("Dynamic_features", backref=backref('lesion', order_by=lesion_id))
    f_morphology = relationship("Morpho_features", backref=backref('lesion', order_by=lesion_id))
    f_texture = relationship("Texture_features", backref=backref('lesion', order_by=lesion_id))
    f_T2 = relationship("T2_features", backref=backref('lesion', order_by=lesion_id))
    annotation = relationship("Annot_record", backref=backref('lesion', order_by=lesion_id))
    segmentation = relationship("Segment_record", backref=backref('lesion', order_by=lesion_id))
        
    def __init__(self, lesionfile, cad_id, dicom_no, accession_no, exam_date, cad_status, mutation, mass_yn, nonmass_yn, finding_side, proc_id, proc_date, proc_side, proc_source, proc_guid, proc_type, lesion_comments, original_report, curve_int, dce_init, dce_delay, label, diagnosis):      
        self.lesionfile = lesionfile
        self.cad_pt_no_txt = cad_id
        self.exam_img_dicom_txt = dicom_no
        self.exam_a_number_txt = accession_no
        self.exam_dt_datetime = exam_date
        self.exam_mri_cad_status_txt = cad_status
        self.cad_latest_mutation_status_int = mutation
        self.exam_find_mri_mass_yn = mass_yn
        self.exam_find_mri_nonmass_yn = nonmass_yn
        self.exam_find_side_int = finding_side
        self.proc_pt_procedure_id = proc_id
        self.proc_proc_dt_datetime = proc_date
        self.proc_proc_side_int = proc_side
        self.proc_proc_source_int = proc_source
        self.proc_proc_guid_int = proc_guid
        self.proc_proc_tp_int = proc_type
        self.proc_lesion_comments_txt = lesion_comments
        self.proc_original_report_txt = original_report
        self.find_curve_int = curve_int
        self.find_mri_dce_init_enh_int = dce_init
        self.find_mri_dce_delay_enh_int = dce_delay
        self.lesion_label = label
        self.lesion_diagnosis = diagnosis
        
    def __repr__(self):
        return "<Lesion_record(lesion_id='%s', cad_pt='%s', dicom='%s', a_number='%s', exam_date='%s')>" % (self.lesion_id, self.cad_pt_no_txt, self.exam_img_dicom_txt, self.exam_a_number_txt, self.exam_dt_datetime)


#  created a Mass_record mapping 
class Mass_record(Base):
    """Base for mass_lesion class using Declarative. for table mass_lesion
    attributes:
        self.lesion_id = lesion_id
        self.BenignNMaligNAnt = BenignNMaligNAnt
        self.DynSeries_id = DynSeries_id
        self.T2Series_id = T2Series_id
        self.find_mammo_n_mri_mass_shape_int = mri_mass_shape
        self.find_mri_mass_margin_int = mri_mass_margin
    """
    __tablename__ = 'mass_lesion'
    __table_args__ = {'autoload':True}
    mass_id = Column(Integer, primary_key=True)
    lesion_id = Column(Integer, ForeignKey('lesion.lesion_id'))
    
    def __init__(self, lesion_id, BenignNMaligNAnt, DynSeries_id, T2Series_id, mri_mass_shape, mri_mass_margin):
        self.lesion_id = lesion_id
        self.BenignNMaligNAnt = BenignNMaligNAnt
        self.DynSeries_id = DynSeries_id
        self.T2Series_id = T2Series_id
        self.find_mammo_n_mri_mass_shape_int = mri_mass_shape
        self.find_mri_mass_margin_int = mri_mass_margin
        
    def __repr__(self):
        return "<Mass_record(lesion_id='%s', BenignNMaligNAnt='%s', DynSeries_id='%s', T2Series_id='%s', mri_mass_shape='%s', , mri_mass_margin='%s')>" % (self.lesion_id, self.BenignNMaligNAnt, self.DynSeries_id, self.T2Series_id, self.find_mammo_n_mri_mass_shape_int, self.find_mri_mass_margin_int)


#  created a Mass_record mapping 
class Nonmass_record(Base):
    """Base for mass_lesion class using Declarative. for table mass_lesion
    attributes:
        self.lesion_id = lesion_id
        self.BenignNMaligNAnt = BenignNMaligNAnt
        self.DynSeries_id = DynSeries_id
        self.T2Series_id = T2Series_id
        self.find_mri_nonmass_dist_int = mri_nonmass_dist
        self.find_mri_nonmass_int_enh_int = mri_nonmass_int_enh   
    """
    __tablename__ = 'nonmass_lesion'
    __table_args__ = {'autoload':True}
    nonmass_id = Column(Integer, primary_key=True)
    lesion_id = Column(Integer, ForeignKey('lesion.lesion_id'))
        
    def __init__(self, lesion_id, BenignNMaligNAnt, DynSeries_id, T2Series_id, mri_nonmass_dist, mri_nonmass_int_enh):
        self.lesion_id = lesion_id
        self.BenignNMaligNAnt = BenignNMaligNAnt
        self.DynSeries_id = DynSeries_id
        self.T2Series_id = T2Series_id
        self.find_mri_nonmass_dist_int = mri_nonmass_dist
        self.find_mri_nonmass_int_enh_int = mri_nonmass_int_enh
        
    def __repr__(self):
        return "<Nonmass_record(lesion_id='%s', BenignNMaligNAnt='%s', DynSeries_id='%s', T2Series_id='%s', mri_nonmass_dist='%s', , mri_nonmass_int_enh='%s')>" % (self.lesion_id, self.BenignNMaligNAnt, self.DynSeries_id, self.T2Series_id, self.find_mri_nonmass_dist_int, self.find_mri_nonmass_int_enh_int)


#  created a Cad_record mapping 
class Dynamic_features(Base):
    """Base for Dynamic_features class using Declarative. for table f_dynamic"""
    __tablename__ = 'f_dynamic'
    __table_args__ = {'autoload':True}
    f_dyn_id = Column(Integer, primary_key=True)
    lesion_id = Column(Integer, ForeignKey('lesion.lesion_id'))
    
        
    def __init__(self, lesion_id, A_inside, alpha_inside, beta_inside, iAUC1_inside, Slope_ini_inside, Tpeak_inside,
                 Kpeak_inside, SER_inside, maxCr_inside, peakCr_inside, UptakeRate_inside, washoutRate_inside, maxVr_inside,
                 peakVr_inside, Vr_increasingRate_inside, Vr_decreasingRate_inside, Vr_post_1_inside,
                 A_countor, alpha_countor, beta_countor, iAUC1_countor, Slope_ini_countor, Tpeak_countor,
                 Kpeak_countor, SER_countor, maxCr_countor, peakCr_countor, UptakeRate_countor, washoutRate_countor, maxVr_countor,
                 peakVr_countor, Vr_increasingRate_countor, Vr_decreasingRate_countor, Vr_post_1_countor):
        self.lesion_id = lesion_id
        self.A_inside = A_inside
        self.alpha_inside = alpha_inside
        self.beta_inside = beta_inside
        self.iAUC1_inside = iAUC1_inside
        self.Slope_ini_inside = Slope_ini_inside
        self.Tpeak_inside = Tpeak_inside
        self.Kpeak_inside = Kpeak_inside
        self.SER_inside = SER_inside
        self.maxCr_inside = maxCr_inside
        self.peakCr_inside = peakCr_inside
        self.UptakeRate_inside = UptakeRate_inside
        self.washoutRate_inside = washoutRate_inside
        self.maxVr_inside = maxVr_inside
        self.peakVr_inside = peakVr_inside
        self.Vr_increasingRate_inside = Vr_increasingRate_inside
        self.Vr_decreasingRate_inside = Vr_decreasingRate_inside
        self.Vr_post_1_inside = Vr_post_1_inside
        self.A_countor = A_countor
        self.alpha_countor = alpha_countor
        self.beta_countor = beta_countor
        self.iAUC1_countor = iAUC1_countor
        self.Slope_ini_countor = Slope_ini_countor
        self.Tpeak_countor = Tpeak_countor
        self.Kpeak_countor = Kpeak_countor
        self.SER_countor = SER_countor
        self.maxCr_countor = maxCr_countor
        self.peakCr_countor = peakCr_countor
        self.UptakeRate_countor = UptakeRate_countor
        self.washoutRate_countor = washoutRate_countor
        self.maxVr_countor = maxVr_countor
        self.peakVr_countor = peakVr_countor
        self.Vr_increasingRate_countor = Vr_increasingRate_countor
        self.Vr_decreasingRate_countor = Vr_decreasingRate_countor
        self.Vr_post_1_countor = Vr_post_1_countor
        
        
#  created a Cad_record mapping 
class Morpho_features(Base):
    """Base for Morpho_features class using Declarative. for table f_morphology"""
    __tablename__ = 'f_morphology'
    __table_args__ = {'autoload':True}
    f_morpho_id = Column(Integer, primary_key=True)
    lesion_id = Column(Integer, ForeignKey('lesion.lesion_id'))
            
    def __init__(self, lesion_id, min_F_r_i, max_F_r_i, mean_F_r_i, var_F_r_i, skew_F_r_i, kurt_F_r_i,
                 iMax_Variance_uptake, iiMin_change_Variance_uptake, iiiMax_Margin_Gradient, k_Max_Margin_Grad,
                 ivVariance, circularity, irregularity, edge_sharp_mean, edge_sharp_std, max_RGH_mean, max_RGH_mean_k, max_RGH_var, max_RGH_var_k):
        self.lesion_id = lesion_id
        self.min_F_r_i = min_F_r_i
        self.max_F_r_i = max_F_r_i
        self.mean_F_r_i = mean_F_r_i
        self.var_F_r_i = var_F_r_i
        self.skew_F_r_i = skew_F_r_i
        self.kurt_F_r_i = kurt_F_r_i
        self.iMax_Variance_uptake = iMax_Variance_uptake
        self.iiMin_change_Variance_uptake = iiMin_change_Variance_uptake
        self.iiiMax_Margin_Gradient = iiiMax_Margin_Gradient
        self.k_Max_Margin_Grad = k_Max_Margin_Grad
        self.ivVariance = ivVariance
        self.circularity = circularity
        self.irregularity = irregularity
        self.edge_sharp_mean = edge_sharp_mean
        self.edge_sharp_std = edge_sharp_std
        self.max_RGH_mean = max_RGH_mean
        self.max_RGH_mean_k = max_RGH_mean_k
        self.max_RGH_var = max_RGH_var
        self.max_RGH_var_k = max_RGH_var_k        
        

#  created a Texture_features mapping 
class Texture_features(Base):
    """Base for Texture_features class using Declarative. for table f_texture"""
    __tablename__ = 'f_texture'
    __table_args__ = {'autoload':True}
    f_texture_id = Column(Integer, primary_key=True)
    lesion_id = Column(Integer, ForeignKey('lesion.lesion_id'))
        
    def __init__(self, lesion_id, texture_contrast_zero, texture_contrast_quarterRad, texture_contrast_halfRad, texture_contrast_threeQuaRad,  
                 texture_homogeneity_zero, texture_homogeneity_quarterRad, texture_homogeneity_halfRad, texture_homogeneity_threeQuaRad,
                 texture_dissimilarity_zero, texture_dissimilarity_quarterRad, texture_dissimilarity_halfRad, texture_dissimilarity_threeQuaRad,
                 texture_correlation_zero, texture_correlation_quarterRad, texture_correlation_halfRad, texture_correlation_threeQuaRad,
                 texture_ASM_zero, texture_ASM_quarterRad, texture_ASM_halfRad, texture_ASM_threeQuaRad,
                 texture_energy_zero, texture_energy_quarterRad, texture_energy_halfRad, texture_energy_threeQuaRad):
        self.lesion_id = lesion_id   
        self.texture_contrast_zero= texture_contrast_zero
        self.texture_contrast_quarterRad= texture_contrast_quarterRad
        self.texture_contrast_halfRad= texture_contrast_halfRad
        self.texture_contrast_threeQuaRad= texture_contrast_threeQuaRad
        self.texture_homogeneity_zero= texture_homogeneity_zero 
        self.texture_homogeneity_quarterRad= texture_homogeneity_quarterRad
        self.texture_homogeneity_halfRad= texture_homogeneity_halfRad
        self.texture_homogeneity_threeQuaRad= texture_homogeneity_threeQuaRad
        self.texture_dissimilarity_zero= texture_dissimilarity_zero
        self.texture_dissimilarity_quarterRad= texture_dissimilarity_quarterRad
        self.texture_dissimilarity_halfRad=texture_dissimilarity_halfRad
        self.texture_dissimilarity_threeQuaRad=texture_dissimilarity_threeQuaRad
        self.texture_correlation_zero= texture_correlation_zero
        self.texture_correlation_quarterRad= texture_correlation_quarterRad
        self.texture_correlation_halfRad= texture_correlation_halfRad
        self.texture_correlation_threeQuaRad= texture_correlation_threeQuaRad
        self.texture_ASM_zero= texture_ASM_zero
        self.texture_ASM_quarterRad= texture_ASM_quarterRad
        self.texture_ASM_halfRad= texture_ASM_halfRad
        self.texture_ASM_threeQuaRad= texture_ASM_threeQuaRad
        self.texture_energy_zero= texture_energy_zero
        self.texture_energy_quarterRad= texture_energy_quarterRad
        self.texture_energy_halfRad= texture_energy_halfRad
        self.texture_energy_threeQuaRad= texture_energy_threeQuaRad
        


#  created a T2_features mapping 
class T2_features(Base):
    """Base for T2_features class using Declarative. for table f_T2"""
    __tablename__ = 'f_T2'
    __table_args__ = {'autoload':True}
    f_T2_id = Column(Integer, primary_key=True)
    lesion_id = Column(Integer, ForeignKey('lesion.lesion_id'))
            
    def __init__(self, lesion_id, find_t2_signal_int, T2dims, T2spacing, T2fatsat,  
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
        self.lesion_id = lesion_id   
        self.find_t2_signal_int= find_t2_signal_int
        self.T2dims= T2dims
        self.T2spacing= T2spacing
        self.T2fatsat= T2fatsat
        self.T2_muscleSI= T2_muscleSI 
        self.T2_muscleSIstd= T2_muscleSIstd
        self.muscle_scalar_range= muscle_scalar_range
        self.bounds_muscleSI= bounds_muscleSI
        self.T2_lesionSI= T2_lesionSI
        self.T2_lesionSIstd= T2_lesionSIstd
        self.lesion_scalar_range= lesion_scalar_range
        self.LMSIR= LMSIR
        self.T2min_F_r_i= T2min_F_r_i
        self.T2max_F_r_i= T2max_F_r_i
        self.T2mean_F_r_i= T2mean_F_r_i
        self.T2var_F_r_i= T2var_F_r_i
        self.T2skew_F_r_i= T2skew_F_r_i
        self.T2kurt_F_r_i= T2kurt_F_r_i
        self.T2grad_margin= T2grad_margin
        self.T2grad_margin_var= T2grad_margin_var
        self.T2RGH_mean= T2RGH_mean
        self.T2RGH_var= T2RGH_var
        self.T2texture_contrast_zero= T2texture_contrast_zero
        self.T2texture_contrast_quarterRad= T2texture_contrast_quarterRad
        self.T2texture_contrast_halfRad= T2texture_contrast_halfRad
        self.T2texture_contrast_threeQuaRad= T2texture_contrast_threeQuaRad
        self.T2texture_homogeneity_zero= T2texture_homogeneity_zero 
        self.T2texture_homogeneity_quarterRad= T2texture_homogeneity_quarterRad
        self.T2texture_homogeneity_halfRad= T2texture_homogeneity_halfRad
        self.T2texture_homogeneity_threeQuaRad= T2texture_homogeneity_threeQuaRad
        self.T2texture_dissimilarity_zero= T2texture_dissimilarity_zero
        self.T2texture_dissimilarity_quarterRad= T2texture_dissimilarity_quarterRad
        self.T2texture_dissimilarity_halfRad=T2texture_dissimilarity_halfRad
        self.T2texture_dissimilarity_threeQuaRad=T2texture_dissimilarity_threeQuaRad
        self.T2texture_correlation_zero= T2texture_correlation_zero
        self.T2texture_correlation_quarterRad= T2texture_correlation_quarterRad
        self.T2texture_correlation_halfRad= T2texture_correlation_halfRad
        self.T2texture_correlation_threeQuaRad= T2texture_correlation_threeQuaRad
        self.T2texture_ASM_zero= T2texture_ASM_zero
        self.T2texture_ASM_quarterRad= T2texture_ASM_quarterRad
        self.T2texture_ASM_halfRad= T2texture_ASM_halfRad
        self.T2texture_ASM_threeQuaRad= T2texture_ASM_threeQuaRad
        self.T2texture_energy_zero= T2texture_energy_zero
        self.T2texture_energy_quarterRad= T2texture_energy_quarterRad
        self.T2texture_energy_halfRad= T2texture_energy_halfRad
        self.T2texture_energy_threeQuaRad= T2texture_energy_threeQuaRad



#  created a Annot_record mapping 
class Annot_record(Base):
    """Base for Annot_record class using Declarative. for table annotation
    attributes:
        self.lesion_id = lesion_id
        self.AccessionNumber = AccessionNumber
        self.SeriesDate = SeriesDate
        self.SeriesNumber = SeriesNumber
        self.SliceLocation = SliceLocation
        self.SeriesDescription = SeriesDescription
        self.PatientID = PatientID
        self.StudyID = StudyID
        self.note = note
        self.xi_coord = xi_coord
        self.yi_coord = yi_coord
        self.xf_coord = xf_coord
        self.yf_coord = yf_coord
        self.SeriesInstanceUID = SeriesInstanceUID
    """
    __tablename__ = 'annotation'
    __table_args__ = {'autoload':True}
    annot_id = Column(Integer, primary_key=True)
    lesion_id = Column(Integer, ForeignKey('lesion.lesion_id'))
        
    def __init__(self, lesion_id, AccessionNumber, SeriesDate, SeriesNumber, SliceLocation, SeriesDescription,
                 PatientID, StudyID, SeriesInstanceUID, note, xi_coord, yi_coord, xf_coord, yf_coord,
                 pi_ijk, pi_2display, pf_ijk, pf_2display, eu_dist_mkers, eu_dist_seg):
        self.lesion_id = lesion_id
        self.AccessionNumber = AccessionNumber
        self.SeriesDate = SeriesDate
        self.SeriesNumber = SeriesNumber
        self.SliceLocation = SliceLocation
        self.SeriesDescription = SeriesDescription
        self.PatientID = PatientID
        self.StudyID = StudyID
        self.SeriesInstanceUID = SeriesInstanceUID
        self.note = note
        self.xi_coord = xi_coord
        self.yi_coord = yi_coord
        self.xf_coord = xf_coord
        self.yf_coord = yf_coord
        self.pi_ijk = pi_ijk
        self.pi_2display = pi_2display
        self.pf_ijk = pf_ijk
        self.pf_2display = pf_2display
        self.eu_dist_mkers = eu_dist_mkers
        self.eu_dist_seg = eu_dist_seg
        
        
    def __repr__(self):
        return "<Annot_record(lesion_id='%s', PatientID='%s', StudyID='%s', note='%s', xi_coord='%s, yi_coord='%s', xf_coord='%s', yf_coord='%s')>" % (self.lesion_id, self.PatientID, self.StudyID, self.note, self.xi_coord, self.yi_coord, self.xf_coord, self.yf_coord)


#  created a Annot_record mapping 
class Segment_record(Base):
    """Base for Segment_record class using Declarative. for table segmentation
    attributes:
        self.lesion_id = lesion_id
        self.segm_xmin = segm_xmin
        self.segm_xmax = segm_xmax
        self.segm_ymin = segm_ymin
        self.segm_ymax = segm_ymax
        self.segm_zmin = segm_zmin
        self.segm_zmax = segm_zmax
        self.no_pts = no_pts
        self.voi_vol = voi_vol
        self.voi_surface = voi_surface
        self.lesion_centroid_world = lesion_centroid_world
        self.lesion_centroid_ijk = lesion_centroid_ijk
    """
    __tablename__ = 'segmentation'
    __table_args__ = {'autoload':True}
    segm_id = Column(Integer, primary_key=True)
    lesion_id = Column(Integer, ForeignKey('lesion.lesion_id'))
        
    def __init__(self, lesion_id, segm_xmin, segm_xmax, segm_ymin, segm_ymax, segm_zmin, segm_zmax,
                 no_pts, voi_vol, voi_surface, VOI_efect_diameter, lesion_centroid_world, lesion_centroid_ijk):
        self.lesion_id = lesion_id
        self.segm_xmin = segm_xmin
        self.segm_xmax = segm_xmax
        self.segm_ymin = segm_ymin
        self.segm_ymax = segm_ymax
        self.segm_zmin = segm_zmin
        self.segm_zmax = segm_zmax
        self.no_pts = no_pts
        self.voi_vol = voi_vol
        self.voi_surface = voi_surface
        self.voi_effect_dia = VOI_efect_diameter
        self.lesion_centroid_world = lesion_centroid_world
        self.lesion_centroid_ijk = lesion_centroid_ijk
        
        
    def __repr__(self):
        return "<Segment_record(lesion_id='%s', voi_vol='%s', voi_surface='%s', lesion_centroid_world='%s', lesion_centroid_ijk='%s')>" % (self.lesion_id, self.voi_vol, self.voi_surface, self.lesion_centroid_world, self.lesion_centroid_ijk)
