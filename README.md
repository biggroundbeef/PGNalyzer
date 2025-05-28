# PGNalyzer
-*- DESCRIPTION -*-

This is a python package dedicated to the calibration and analysis of Pandonia Global Network (PGN) data. 

Quick Setup Guide:
1) Download PGNalyzer from GitHub.
2) Unzip all files from zipped folder.
3) Locate "pgnalyzer" folder 
4) Copy/Move pgnalyzer folder into python environment.*

* In today's world you may either be simply running python on your computer and mixing all libraries and python versions you've ever used, 
or you will be in a python environment where your packages are easily accessible. If the former pertains to you, I'd suggest learning how
to implement a python envrionment on your computer. For the latter, simply copy/move the "pgnalyzer" folder into your packages directory.

Example Code: 
from pgnalyzer.io.read_data import Reader
    ron = Reader()
    in_PGN_data_test = path
    in_model_path = path
    routines, dora_data = ron.read_PGN_txt()
    corrected_counts = ron.get_dark_corrected_counts()
    spectral_calibrated_channels, a, b, corr = ron.get_spectral_calibrated_channels()
    ron.plot_corrected_data()

Current Capabilites:
    Objects:
    -> Reader()
    This is found in pgnalyzer.io.read_data
        -> Reader.read_PGN_txt() 
        -> Reader.get_dark_corrected_counts()
        -> Reader.get_spectral_calibrated_channels()
        -> Reader.plot_corrected_data()

    Modules:
    astro
    -> sza 
        -> get_solar_zenith_angle()
    
    dark_funcs
    -> raw_text
        -> is_dark()
        -> simplify_darks()
        -> find_best_dark() 

    io
    -> model 
        -> get_model()
    -> read_data
        *** This file currently houses the Reader class. Reader class will be properly created in future updates.

    mods 
    -> data 
        -> smooth()
    -> dates 
        -> update_raw()
        -> update_str()
        -> ct2lst
    -> raw_text
        -> find_data_start()
        -> find_intro_header()
        -> find_lat_long()
        -> find_headers()
        -> sort_routine()
        -> separate_routines()
        -> read_PGN()
        -> apply_dark_correction()

    spec_cal 
    -> solve
        -> get_channel() 
        -> fit_spectrum()


PGNalyzer 
Version 1.0 alpha
Release Date: 05/28/2025

Creator:
Eduardo Sanchez Jr. (New Mexico State University)

Supporters:
Dr. Juie Shetye
Dr. Dong Wu
Peter Pantina
Lipi Mukherjee
