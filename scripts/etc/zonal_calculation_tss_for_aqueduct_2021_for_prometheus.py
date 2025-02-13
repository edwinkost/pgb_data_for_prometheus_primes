#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import datetime
import calendar
import time
import re
import subprocess
import pcraster as pcr
import netCDF4 as nc
import numpy as np
import virtualOS as vos

class MakingNetCDF():
    
    def __init__(self, cloneMapFile, attribute = None, cellSizeInArcMinutes = None):
        		
        # cloneMap
        # - the cloneMap must be at 5 arc min resolution
        cloneMap = pcr.readmap(cloneMapFile)
        cloneMap = pcr.boolean(1.0)
        
        # latitudes and longitudes
        self.latitudes  = np.unique(pcr.pcr2numpy(pcr.ycoordinate(cloneMap), vos.MV))[::-1]
        self.longitudes = np.unique(pcr.pcr2numpy(pcr.xcoordinate(cloneMap), vos.MV))

        #~ # properties of the clone map
        #~ # - number of rows and columns
        #~ self.nrRows       = np.round(pcr.clone().nrRows())    
        #~ self.nrCols       = np.round(pcr.clone().nrCols())  
        #~ # - upper right coordinate, unit: arc degree ; must be integer (without decimals)
        #~ self.minLongitude = np.round(pcr.clone().west() , 0)         
        #~ self.maxLatitude  = np.round(pcr.clone().north(), 0)
        #~ # - cell resolution, unit: arc degree
        #~ self.cellSize     = pcr.clone().cellSize()
        #~ if cellSizeInArcMinutes != None: self.cellSize = cellSizeInArcMinutes / 60.0 
        #~ # - lower right coordinate, unit: arc degree ; must be integer (without decimals)
        #~ self.maxLongitude = np.round(self.minLongitude + self.cellSize*self.nrCols, 0)         
        #~ self.minLatitude  = np.round(self.maxLatitude  - self.cellSize*self.nrRows, 0)
        #~ 
        #~ # latitudes and longitudes for netcdf files
        #~ latMin = self.minLatitude  + self.cellSize / 2
        #~ latMax = self.maxLatitude  - self.cellSize / 2
        #~ lonMin = self.minLongitude + self.cellSize / 2
        #~ lonMax = self.maxLongitude - self.cellSize / 2
        #~ self.longitudes = np.arange(lonMin,lonMax+self.cellSize, self.cellSize)
        #~ self.latitudes=   np.arange(latMax,latMin-self.cellSize,-self.cellSize)
        
        # netCDF format and attributes:
        self.format = 'NETCDF4'
        self.attributeDictionary = {}
        if attribute == None:
            self.attributeDictionary['institution'] = "None"
            self.attributeDictionary['title'      ] = "None"
            self.attributeDictionary['description'] = "None"
        else:
            self.attributeDictionary = attribute

    def createNetCDF(self,ncFileName,varName,varUnit):

        rootgrp= nc.Dataset(ncFileName,'w',format= self.format)

        #-create dimensions - time is unlimited, others are fixed
        rootgrp.createDimension('time',None)
        rootgrp.createDimension('lat',len(self.latitudes))
        rootgrp.createDimension('lon',len(self.longitudes))

        date_time= rootgrp.createVariable('time','f4',('time',))
        date_time.standard_name= 'time'
        date_time.long_name= 'Days since 1901-01-01'

        date_time.units= 'Days since 1901-01-01' 
        date_time.calendar= 'standard'

        lat= rootgrp.createVariable('lat','f4',('lat',))
        lat.long_name= 'latitude'
        lat.units= 'degrees_north'
        lat.standard_name = 'latitude'

        lon= rootgrp.createVariable('lon','f4',('lon',))
        lon.standard_name= 'longitude'
        lon.long_name= 'longitude'
        lon.units= 'degrees_east'

        lat[:]= self.latitudes
        lon[:]= self.longitudes

        shortVarName = varName
        var= rootgrp.createVariable(shortVarName,'f4',('time','lat','lon',) ,fill_value=vos.MV,zlib=True)
        var.standard_name = shortVarName
        var.long_name = shortVarName
        var.units = varUnit

        attributeDictionary = self.attributeDictionary
        for k, v in attributeDictionary.items():
          setattr(rootgrp,k,v)

        rootgrp.sync()
        rootgrp.close()

    def writePCR2NetCDF(self,ncFileName,varName,varField,timeStamp,posCnt):

        #-write data to netCDF
        rootgrp= nc.Dataset(ncFileName,'a')    

        shortVarName= varName        

        date_time= rootgrp.variables['time']
        date_time[posCnt]= nc.date2num(timeStamp,date_time.units,date_time.calendar)

        rootgrp.variables[shortVarName][posCnt,:,:]= (varField)

        rootgrp.sync()
        rootgrp.close()

if __name__ == "__main__":
    
    # clone, landmask and cell area files
    landmask05minFile    = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/lddsound_05min_version_20210330.map"
    #~ landmask05minFile = "/projects/0/dfguu/data/hydroworld/others/RhineMeuse/RhineMeuse05min.landmask.map"
    cloneMapFileName     = landmask05minFile 
    cellSizeInArcMinutes = 5.0 
    cellArea05minFile    = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/cdo_gridarea_clone_global_05min_correct_lats.nc"
    # set clone
    pcr.setclone(landmask05minFile)

    # ~ edwinbar@fcn73.local.snellius.surf.nl:/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/historical_and_ssp_files$ ls -lah
    # ~ total 614G
    # ~ drwxr-xr-x 2 edwin edwin 4.0K Jan 21 16:19 .
    # ~ drwxr-xr-x 5 edwin edwin 4.0K Apr  1  2022 ..
    # ~ -rw-r--r-- 1 edwin edwin  51G Apr  1  2022 domestic_water_demand_historical_1960-2019.nc
    # ~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 domestic_water_demand_ssp1_2000-2100.nc
    # ~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 domestic_water_demand_ssp3_2000-2100.nc
    # ~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 domestic_water_demand_ssp5_2000-2100.nc
    # ~ -rw-r--r-- 1 edwin edwin  51G Apr  1  2022 industry_water_demand_historical_1960-2019.nc
    # ~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 industry_water_demand_ssp1_2000-2100.nc
    # ~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 industry_water_demand_ssp3_2000-2100.nc
    # ~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 industry_water_demand_ssp5_2000-2100.nc

   
    # class (country) ids
    uniqueIDsFile = "/home/edwinbar/github/edwinkost/pgb_data_for_prometheus_primes/iso3_countries/row_number_CNTR_RG_01M_2020_4326.shp.tif.map"

    # output directory of this analusis
    outputDirectory = "/scratch-shared/edwinbar/pgb_for_prometheus/test/"
    outputDirectory = sys.argv[1]
    
    # output file code/pattern
    output_file_code = "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference"
    output_file_code = sys.argv[2] 
    
    # codes for gcm, 'scenario', etc
    gcmCode      = sys.argv[2]
    scenarioCode = sys.argv[3]
    output_file_code = "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s" % (gcmCode, scenarioCode) 

    # start year and end year
    staYear = 2000
    endYear = 2005
    staYear = int(sys.argv[4])
    endYear = int(sys.argv[5])

    # input files
    #
    wt
    inputFiles['industryGrossDemand'] = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_totalRunoff_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))

    # - output directory of PCR-GLOBWB run:
    # ~ pgb_output_dir  = "/depfg/sutan101/pcrglobwb_wri_aqueduct_2021/pcrglobwb_aqueduct_2021_monthly_annual_files/version_2021-09-16_merged/gswp3-w5e5/historical-reference/"
    pgb_output_dir  = "/depfg/sutan101/pcrglobwb_wri_aqueduct_2021/pcrglobwb_aqueduct_2021_monthly_annual_files/version_2021-09-16_merged/%s/%s/" % (gcmCode, scenarioCode)
    #
    inputFiles = {}
    #
    # - unit of the following is m.year and flux values given are over the entire cell area
    #
    # ~ inputFiles['totalRunoff']                  = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_totalRunoff_global_yearly-total_1960_2019_basetier1.nc"
    # ~ inputFiles['gwRecharge']                   = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_gwRecharge_global_yearly-total_1960_2019_basetier1.nc"
    # ~ inputFiles['fossilGroundwaterAbstraction'] = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_fossilGroundwaterAbstraction_global_yearly-total_1960_2019_basetier1.nc"
    # ~ inputFiles['desalinationAbstraction']      = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_desalinationAbstraction_global_yearly-total_1960_2019_basetier1.nc"
    # ~ inputFiles['domesticWaterWithdrawal']      = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_domesticWaterWithdrawal_global_yearly-total_1960_2019_basetier1.nc"
    # ~ inputFiles['industryWaterWithdrawal']      = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_industryWaterWithdrawal_global_yearly-total_1960_2019_basetier1.nc"
    # ~ inputFiles['livestockWaterWithdrawal']     = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_livestockWaterWithdrawal_global_yearly-total_1960_2019_basetier1.nc"
    # ~ inputFiles['irrigationWaterWithdrawal']    = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_irrigationWaterWithdrawal_global_yearly-total_1960_2019_basetier1.nc"

    inputFiles['totalRunoff']                  = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_totalRunoff_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))
    inputFiles['gwRecharge']                   = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_gwRecharge_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))
    inputFiles['fossilGroundwaterAbstraction'] = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_fossilGroundwaterAbstraction_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))
    inputFiles['desalinationAbstraction']      = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_desalinationAbstraction_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))
    inputFiles['domesticWaterWithdrawal']      = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_domesticWaterWithdrawal_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))
    inputFiles['industryWaterWithdrawal']      = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_industryWaterWithdrawal_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))
    inputFiles['livestockWaterWithdrawal']     = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_livestockWaterWithdrawal_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))
    inputFiles['irrigationWaterWithdrawal']    = pgb_output_dir + "/" + "pcrglobwb_cmip6-isimip3-%s_image-aqueduct_%s_irrigationWaterWithdrawal_global_yearly-total_%s_%s_basetier1.nc" % (gcmCode, scenarioCode, str(staYear), str(endYear))


    # capri also needs irrigation efficiency
    inputFiles['irrigationEfficiency']         = irrigation_eff_file 
    
    
    # output that will be calculated 
    output = {}
    variable_names  = inputFiles.keys()
    for var in variable_names:
        output[var] = {}
        # ~ output[var]['file_name'] = outputDirectory + "/" + output_file_code + "_" + str(var) + "_annual_country.nc"
        output[var]['file_name'] = outputDirectory + "/" + output_file_code + "_" + str(var) + "_annual_country.nc"
        output[var]['unit']      = "km3.year-1"
        if var == "irrigationEfficiency": output[var]['unit'] = "1" 
        output[var]['pcr_value'] = None

    # making output and temporary directories
    if os.path.exists(outputDirectory):
        shutil.rmtree(outputDirectory)
    os.makedirs(outputDirectory)
    # - moving to the output directory
    os.chdir(outputDirectory)
    # - temporary directory
    tmp_directory = outputDirectory + "/tmp/"
    os.makedirs(tmp_directory)
    # - table directory
    table_directory = outputDirectory + "/table/"
    os.makedirs(table_directory)
    
    # attribute for netCDF files 
    attributeDictionary = {}
    attributeDictionary['title'      ]  = "PCR-GLOBWB output files"
    attributeDictionary['institution']  = "Dept. of Physical Geography, Utrecht University"
    attributeDictionary['source'     ]  = pgb_output_dir
    attributeDictionary['history'    ]  = "None"
    attributeDictionary['references' ]  = "None"
    attributeDictionary['comment'    ]  = "The conuntry code is based on the " +  str(uniqueIDsFile) + "."
    # additional attribute defined in PCR-GLOBWB 
    attributeDictionary['description']  = "PCR-GLOBWB output files from the Aqueduct run 2021, prepared by Edwin H. Sutanudjaja."

    # initiate the netcd object: 
    tssNetCDF = MakingNetCDF(cloneMapFile = cloneMapFileName, \
                             attribute = attributeDictionary, \
                             cellSizeInArcMinutes = cellSizeInArcMinutes)
    # making netcdf files:
    for var in variable_names:
        tssNetCDF.createNetCDF(output[var]['file_name'], var, output[var]['unit'])

    # class (country) ids
    uniqueIDs = pcr.nominal(\
                vos.readPCRmapClone(uniqueIDsFile, cloneMapFileName, tmp_directory, 
                                    None, False, None, True))
    uniqueIDs = pcr.ifthen(pcr.scalar(uniqueIDs) >= 0.0, uniqueIDs)
    
    # landmask                               
    landmask = pcr.defined(pcr.readmap(landmask05minFile))
    landmask = pcr.ifthen(landmask, landmask)
    # - extending landmask with uniqueIDs
    landmask = pcr.cover(landmask, pcr.defined(uniqueIDs))
    
    # extending class (country) ids
    max_step = 7
    for i in range(1, max_step+1, 1):
        cmd = "Extending class: step "+str(i)+" from " + str(max_step)
        print(cmd)
        uniqueIDs = pcr.cover(uniqueIDs, pcr.windowmajority(uniqueIDs, 0.5))
    # - cover the rest with a new id
    uniqueIDs = pcr.cover(uniqueIDs, pcr.nominal(pcr.mapmaximum(pcr.scalar(uniqueIDs)) + 1000))
    # - use only cells within the landmask
    uniqueIDs = pcr.ifthen(landmask, uniqueIDs)
    pcr.report(uniqueIDs, "class_ids.map")                                
    
    # cell area at 5 arc min resolution
    cellArea = vos.readPCRmapClone(cellArea05minFile,
                                   cloneMapFileName, tmp_directory)
    cellArea = pcr.ifthen(landmask, cellArea)
    
    # get a sample cell for every id
    x_min_for_each_id = pcr.areaminimum(pcr.xcoordinate(pcr.boolean(1.0)), uniqueIDs)
    sample_cells      = pcr.xcoordinate(pcr.boolean(1.0)) == x_min_for_each_id
    y_min_for_each_id = pcr.areaminimum(pcr.ycoordinate(sample_cells), uniqueIDs)
    sample_cells      = pcr.ycoordinate(sample_cells) == y_min_for_each_id
    uniqueIDs_sample  = pcr.ifthen(sample_cells, uniqueIDs)
    # - save it to a pcraster map file
    pcr.report(uniqueIDs_sample, "sample.ids")                                

    # calculate the country values 
    index = 0 # for posCnt
    for iYear in range(staYear,endYear+1):
        
        # time stamp and index for netcdf files:
        index = index + 1
        timeStamp = datetime.datetime(int(iYear), int(12), int(31), int(0))
        fulldate = '%4i-%02i-%02i'  %(int(iYear), int(12), int(31))
        print(fulldate)

        # reading pcraster or netcdf files:
        for var in inputFiles.keys():        
            
            print(inputFiles[var])   

            if var == "irrigationEfficiency":
            
                if iYear == staYear:
                
                    irrigationEfficiency = vos.readPCRmapClone(inputFiles[var],
                                                               cloneMapFileName, tmp_directory)

                    irrigationEfficiency = pcr.cover(irrigationEfficiency, 1.0)
                    irrigationEfficiency = pcr.max(0.1, irrigationEfficiency)

                output[var]['pcr_value'] = irrigationEfficiency


            if var != "irrigationEfficiency":

                # reading PCR-GLOBWB values from netcdf files with time component

                fulldate_for_reading_netcdf = fulldate
                output[var]['pcr_value'] = vos.netcdf2PCRobjClone(ncFile = inputFiles[var],\
                                                                  varName = "automatic",\
                                                                  dateInput = fulldate_for_reading_netcdf,
                                                                  useDoy = None,
                                                                  cloneMapFileName  = cloneMapFileName,
                                                                  LatitudeLongitude = True,
                                                                  specificFillValue = None)
        
        # upscaling to the class (country) units and writing to netcdf files and a table
        for var in output.keys():
            
            print(var)
            
            if var == "irrigationEfficiency":

                # upscaling to the class (country) unit - using areaaveraga
                pcrValue = pcr.areaaverage(output[var]['pcr_value'], uniqueIDs)

            if var != "irrigationEfficiency":

                # covering the map with zero
                pcrValue = pcr.cover(output[var]['pcr_value'], 0.0)
                
                # convert values from m to m3
                pcrValue =  pcrValue * cellArea
			    
                # upscaling to the class (country) units and converting the units to km3/year - using areatotal
                pcrValue = pcr.areatotal(pcrValue, uniqueIDs) / (1000. * 1000. * 1000.)
            
            # write values to a netcdf file
            ncFileName = output[var]['file_name']
            varField = pcr.pcr2numpy(pcrValue, vos.MV)
            tssNetCDF.writePCR2NetCDF(ncFileName, var, varField, timeStamp, posCnt = index - 1)
            
            # plot the values at sample cells only and write values to a temporary pcraster map
            pcrFileName = str(tmp_directory) + "/" + str(var) + ".tmp"
            pcr.report(pcr.ifthen(pcr.defined(uniqueIDs_sample), pcrValue), pcrFileName)

        # write class values to a table
        # - command line to call map2col
        cmd    = 'map2col -x 1 -y 2 -m NA sample.ids'
        # - header for the table
        header = "x y class_id"
        # - txt file that contains the table
        txt_file = open(table_directory + "/" + output_file_code + "_summary_" + str(iYear) + ".txt", "w")
        for var in output.keys():
            header += " " + str(var)
            if var != "irrigationEfficiency": header += "_km3"
            cmd    += " " + str(tmp_directory) + "/" + str(var) + ".tmp"
        cmd += " " + str(tmp_directory) + "/" + output_file_code + "_summary_" + str(iYear) + ".txt.tmp"
        print(cmd)
        os.system(cmd)
        # - add header to txt file
        header += "\n" 
        txt_file.write(header)
        # - add map2col output to the txt_file
        map2col_file = open(tmp_directory+"/" + output_file_code + "_summary_" + str(iYear) + ".txt.tmp", "r")
        txt_file.write(map2col_file.read())
        # - close all open txt files
        txt_file.close()
        map2col_file.close()
        
        # remove all temporary files
        cmd = 'rm -r '+ tmp_directory + "/*"
        print(cmd)
        os.system(cmd)
        

