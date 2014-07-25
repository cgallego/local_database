# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 10:20:41 2014

@author: Cristina Gallego
"""
import sys, os
import string
import datetime
import numpy as np

from sqlalchemy import *
from sqlalchemy.orm import *
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy import MetaData, Table, Column, DateTime, Integer, String, Boolean, Float
from sqlalchemy.orm import mapper, sessionmaker, polymorphic_union

#import mydatabase
from mybase import Base, myengine

# he Declarative base class also contains a catalog of all the Table objects that have been defined called MetaData
metadata = MetaData()

lesion = Table('lesion', metadata,
    Column('lesion_id', Integer, primary_key=True),
    Column('lesionfile', String(50)),
    Column('cad_pt_no_txt', String(50)),
    Column('exam_img_dicom_txt', String(50)),
    Column('exam_a_number_txt', String(50)),
    Column('exam_dt_datetime', DateTime),
    Column('exam_mri_cad_status_txt', String(50)),
    Column('cad_latest_mutation_status_int', String(50)),
    Column('exam_find_mri_mass_yn', Boolean),
    Column('exam_find_mri_nonmass_yn', Boolean),
    Column('exam_find_side_int', String(50)),
    Column('proc_pt_procedure_id', String(50)),
    Column('proc_proc_dt_datetime', DateTime),
    Column('proc_proc_side_int', String(50)),
    Column('proc_proc_source_int', String(50)),
    Column('proc_proc_guid_int', String(50)),
    Column('proc_proc_tp_int', String),
    Column('proc_lesion_comments_txt', String),
    Column('proc_original_report_txt', String),
    Column('find_curve_int', String(50)),
    Column('find_mri_dce_init_enh_int', String(50)),
    Column('find_mri_dce_delay_enh_int', String(50)),
    Column('lesion_label', String(50)),
    Column('lesion_diagnosis', String)    
)

mass_lesion = Table('mass_lesion', metadata,
    Column('mass_id', Integer, primary_key=True),
    Column('lesion_id', ForeignKey('lesion.lesion_id', ondelete="CASCADE")),
    Column('BenignNMaligNAnt', String(40)),
    Column('DynSeries_id', String(40)),
    Column('T2Series_id', String(40)),
    Column('find_mammo_n_mri_mass_shape_int', String(50)),
    Column('find_mri_mass_margin_int', String(50))
)

nonmass_lesion = Table('nonmass_lesion', metadata,
    Column('nonmass_id', Integer, primary_key=True),
    Column('lesion_id', ForeignKey('lesion.lesion_id', ondelete="CASCADE")),
    Column('BenignNMaligNAnt', String(40)),
    Column('DynSeries_id', String(40)),
    Column('T2Series_id', String(40)),
    Column('find_mri_nonmass_dist_int', String(50)),
    Column('find_mri_nonmass_int_enh_int', String(50))
)

f_dynamic = Table('f_dynamic', metadata,
    Column('f_dyn_id', Integer, primary_key=True),
    Column('lesion_id', ForeignKey('lesion.lesion_id', ondelete="CASCADE")),
    Column('A_inside',  Float),
    Column('alpha_inside',  Float),
    Column('beta_inside',  Float),
    Column('iAUC1_inside',  Float),
    Column('Slope_ini_inside',  Float),
    Column('Tpeak_inside',  Float),
    Column('Kpeak_inside',  Float),
    Column('SER_inside',  Float),
    Column('maxCr_inside',  Float),
    Column('peakCr_inside',  Float),
    Column('UptakeRate_inside',  Float),
    Column('washoutRate_inside',  Float),
    Column('maxVr_inside',  Float),
    Column('peakVr_inside',  Float),
    Column('Vr_increasingRate_inside',  Float),
    Column('Vr_decreasingRate_inside',  Float),
    Column('Vr_post_1_inside',  Float),
    Column('A_countor',  Float),
    Column('alpha_countor',  Float),
    Column('beta_countor',  Float),
    Column('iAUC1_countor',  Float),
    Column('Slope_ini_countor',  Float),
    Column('Tpeak_countor',  Float),
    Column('Kpeak_countor',  Float),
    Column('SER_countor',  Float),
    Column('maxCr_countor',  Float),
    Column('peakCr_countor',  Float),
    Column('UptakeRate_countor',  Float),
    Column('washoutRate_countor',  Float),
    Column('maxVr_countor',  Float),
    Column('peakVr_countor',  Float),
    Column('Vr_increasingRate_countor',  Float),
    Column('Vr_decreasingRate_countor',  Float),
    Column('Vr_post_1_countor',  Float)
)

f_morphology = Table('f_morphology', metadata,
    Column('f_morpho_id', Integer, primary_key=True),
    Column('lesion_id', ForeignKey('lesion.lesion_id', ondelete="CASCADE")),
    Column('min_F_r_i', Numeric),
    Column('max_F_r_i', Numeric),
    Column('mean_F_r_i', Numeric),
    Column('var_F_r_i', Numeric),
    Column('skew_F_r_i', Numeric),
    Column('kurt_F_r_i', Numeric),
    Column('iMax_Variance_uptake', Numeric),
    Column('iiMin_change_Variance_uptake', Numeric),
    Column('iiiMax_Margin_Gradient', Numeric),  
    Column('k_Max_Margin_Grad', Numeric),
    Column('ivVariance', Numeric),
    Column('circularity', Numeric),
    Column('irregularity', Numeric),
    Column('edge_sharp_mean', Numeric),
    Column('edge_sharp_std', Numeric),
    Column('max_RGH_mean', Numeric),
    Column('max_RGH_mean_k', Numeric),
    Column('max_RGH_var', Numeric),
    Column('max_RGH_var_k', Numeric)
)

f_texture = Table('f_texture', metadata,
    Column('f_texture_id', Integer, primary_key=True),
    Column('lesion_id', ForeignKey('lesion.lesion_id', ondelete="CASCADE")),
    Column('texture_contrast_zero', Numeric),
    Column('texture_contrast_quarterRad', Numeric),   
    Column('texture_contrast_halfRad', Numeric),
    Column('texture_contrast_threeQuaRad', Numeric),  
    Column('texture_homogeneity_zero', Numeric), 
    Column('texture_homogeneity_quarterRad', Numeric),
    Column('texture_homogeneity_halfRad', Numeric),
    Column('texture_homogeneity_threeQuaRad', Numeric),
    Column('texture_dissimilarity_zero', Numeric),
    Column('texture_dissimilarity_quarterRad', Numeric),
    Column('texture_dissimilarity_halfRad', Numeric),
    Column('texture_dissimilarity_threeQuaRad', Numeric),
    Column('texture_correlation_zero', Numeric),
    Column('texture_correlation_quarterRad', Numeric),
    Column('texture_correlation_halfRad', Numeric),
    Column('texture_correlation_threeQuaRad', Numeric),
    Column('texture_ASM_zero', Numeric),
    Column('texture_ASM_quarterRad', Numeric),
    Column('texture_ASM_halfRad', Numeric),
    Column('texture_ASM_threeQuaRad', Numeric),
    Column('texture_energy_zero', Numeric),
    Column('texture_energy_quarterRad', Numeric),
    Column('texture_energy_halfRad', Numeric),
    Column('texture_energy_threeQuaRad', Numeric)
)

f_T2 = Table('f_T2', metadata,
    Column('f_T2_id', Integer, primary_key=True),
    Column('lesion_id', ForeignKey('lesion.lesion_id', ondelete="CASCADE")),
    Column('find_t2_signal_int', String),
    Column('T2dims', String),
    Column('T2spacing', String),
    Column('T2fatsat', String),
    Column('T2_muscleSI', Numeric),
    Column('T2_muscleSIstd', Numeric),
    Column('muscle_scalar_range', String),
    Column('bounds_muscleSI', String),
    Column('T2_lesionSI', Numeric),
    Column('T2_lesionSIstd', Numeric),
    Column('lesion_scalar_range', String),
    Column('LMSIR', Numeric),
    Column('T2min_F_r_i', Numeric),
    Column('T2max_F_r_i', Numeric),
    Column('T2mean_F_r_i', Numeric),
    Column('T2var_F_r_i', Numeric),
    Column('T2skew_F_r_i', Numeric),
    Column('T2kurt_F_r_i', Numeric),
    Column('T2grad_margin', Numeric),
    Column('T2grad_margin_var', Numeric),
    Column('T2RGH_mean', Numeric),
    Column('T2RGH_var', Numeric),
    Column('T2texture_contrast_zero', Numeric),
    Column('T2texture_contrast_quarterRad', Numeric),   
    Column('T2texture_contrast_halfRad', Numeric),
    Column('T2texture_contrast_threeQuaRad', Numeric),  
    Column('T2texture_homogeneity_zero', Numeric), 
    Column('T2texture_homogeneity_quarterRad', Numeric),
    Column('T2texture_homogeneity_halfRad', Numeric),
    Column('T2texture_homogeneity_threeQuaRad', Numeric),
    Column('T2texture_dissimilarity_zero', Numeric),
    Column('T2texture_dissimilarity_quarterRad', Numeric),
    Column('T2texture_dissimilarity_halfRad', Numeric),
    Column('T2texture_dissimilarity_threeQuaRad', Numeric),
    Column('T2texture_correlation_zero', Numeric),
    Column('T2texture_correlation_quarterRad', Numeric),
    Column('T2texture_correlation_halfRad', Numeric),
    Column('T2texture_correlation_threeQuaRad', Numeric),
    Column('T2texture_ASM_zero', Numeric),
    Column('T2texture_ASM_quarterRad', Numeric),
    Column('T2texture_ASM_halfRad', Numeric),
    Column('T2texture_ASM_threeQuaRad', Numeric),
    Column('T2texture_energy_zero', Numeric),
    Column('T2texture_energy_quarterRad', Numeric),
    Column('T2texture_energy_halfRad', Numeric),
    Column('T2texture_energy_threeQuaRad', Numeric)
)
            

annotation = Table('annotation', metadata,
    Column('annot_id', Integer, primary_key=True),
    Column('lesion_id', ForeignKey('lesion.lesion_id', ondelete="CASCADE")),
    Column('AccessionNumber', String(40)),
    Column('SeriesDate', String(40)),
    Column('SeriesNumber', String(40)),
    Column('SliceLocation', String(50)),
    Column('SeriesDescription', String(50)),
    Column('PatientID', String(50)),
    Column('DicomExamNumber', String(50)),
    Column('SeriesInstanceUID', String),
    Column('note', String(50)),
    Column('xi_coord', Integer), 
    Column('yi_coord', Integer),
    Column('xf_coord', Integer), 
    Column('yf_coord', Integer),
    Column('pi_ijk', String(50)),
    Column('pi_2display', String(50)),
    Column('pf_ijk', String(50)),
    Column('pf_2display', String(50)),
    Column('eu_dist_mkers', Integer),
    Column('eu_dist_seg', Integer)
)

annotation = Table('segmentation', metadata,
    Column('segm_id', Integer, primary_key=True),
    Column('lesion_id', ForeignKey('lesion.lesion_id', ondelete="CASCADE")),
    Column('segm_xmin', Numeric),
    Column('segm_xmax', Numeric),
    Column('segm_ymin', Numeric),
    Column('segm_ymax', Numeric),   
    Column('segm_zmin', Numeric),
    Column('segm_zmax', Numeric),
    Column('no_pts', Numeric),
    Column('voi_vol', Numeric),
    Column('voi_surface', Numeric),
    Column('voi_effect_dia', Numeric),
    Column('lesion_centroid_world', String(50)),
    Column('lesion_centroid_ijk', String(50))   
)


        
# configure myengine and create tables with desired options
metadata.create_all(myengine)



