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
    
    def __init__(self, cloneMapFile, attribute=None, cellSizeInArcMinutes=None):
        		
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
    landmask05minFile    = "/projects/0/dfguu/data/hydroworld/PCRGLOBWB20/input5min/routing/lddsound_05min.map"
    #~ landmask05minFile = "/projects/0/dfguu/data/hydroworld/others/RhineMeuse/RhineMeuse05min.landmask.map"
    cloneMapFileName     = landmask05minFile 
    cellSizeInArcMinutes = 5.0 
    cellArea05minFile    = "/projects/0/dfguu/data/hydroworld/PCRGLOBWB20/input5min/routing/cellsize05min.correct.map"
    # set clone
    pcr.setclone(landmask05minFile)
    
    # output directory
    outputDirectory = "/scratch-shared/edwinhs/country_water_use_and_demand_for_pcrglobwb2.0_paper/"
    
    # start year and end year
    staYear = 1958
    endYear = 2015

    # input files
    #
    # - output directory of PCR-GLOBWB 2:
    inputDirectory  = "/scratch-shared/edwin/05min_runs_for_land_subsidence_deltares_march_april_2017/4LCs_accutraveltime_cru-forcing/non-natural/merged_1958_to_2015/global/netcdf/"
    # - water demand directory (domestic, industrial, livestock; no aggricultural)
    waterDemandDirectory = "/scratch-shared/edwinhs/05min_runs_for_gmd_paper_30_oct_2017/05min_runs_4LCs_accutraveltime_cru-forcing_1958-2015/non-natural_starting_from_1958/analysis/water_demand/annual_average_m_per_day/" 
    #
    inputFiles = {}
    #
    # - unit of the following is m.year and flux values given are over the entire cell area
    #
    inputFiles["domesticGrossDemand"   ]      = waterDemandDirectory + "/" + "domesticGrossDemand.nc"
    inputFiles["industryGrossDemand"   ]      = waterDemandDirectory + "/" + "industryGrossDemand.nc"
    inputFiles["livestockGrossDemand"  ]      = waterDemandDirectory + "/" + "livestockGrossDemand.nc"
    #                                                                        
    inputFiles["domesticNettoDemand"   ]      = waterDemandDirectory + "/" + "domesticNettoDemand.nc"
    inputFiles["industryNettoDemand"   ]      = waterDemandDirectory + "/" + "industryNettoDemand.nc"
    inputFiles["livestockNettoDemand"  ]      = waterDemandDirectory + "/" + "livestockNettoDemand.nc"
    #
    #~ inputFiles["domestic_water_withdrawal"     ]    = inputDirectory + "/" + "domesticWaterWithdrawal_annuaTot_output"
    #~ inputFiles["industry_water_withdrawal"     ]    = inputDirectory + "/" + "industryWaterWithdrawal_annuaTot_output"
    #~ inputFiles["livestock_water_withdrawal"    ]    = inputDirectory + "/" + "livestockWaterWithdrawal_annuaTot_output"
    #~ #
    #~ inputFiles["non_irrigation_consumption"    ]    = inputDirectory + "/" + "nonIrrWaterConsumption_annuaTot_output"
    #~ #
    #~ inputFiles["precipitation"                 ]    = inputDirectory + "/" + "precipitation_annuaTot_output"
    #~ inputFiles["total_runoff"                  ]    = inputDirectory + "/" + "totalRunoff_annuaTot_output"
    #~ inputFiles["total_evaporation"             ]    = inputDirectory + "/" + "totalEvaporation_annuaTot_output"
    #~ inputFiles["total_groundwater_recharge"    ]    = inputDirectory + "/" + "gwRecharge_annuaTot_output"
    #~ #
    #~ inputFiles["total_abstraction"             ]    = inputDirectory + "/" + "totalAbstraction_annuaTot_output"
    #~ inputFiles["desalination_abstraction"      ]    = inputDirectory + "/" + "desalinationAbstraction_annuaTot_output"
    #~ inputFiles["surface_water_abstraction"     ]    = inputDirectory + "/" + "surfaceWaterAbstraction_annuaTot_output"
    #~ inputFiles["total_groundwater_abstraction" ]    = inputDirectory + "/" + "totalGroundwaterAbstraction_annuaTot_output_1958-12-31_to_2015-12-31.nc"
    #~ inputFiles["fossil_groundwater_abstraction"]    = inputDirectory + "/" + "fossilGroundwaterAbstraction_annuaTot_output"
    #
    #~ # - unit of the following is m.year and flux values given are over the entire cell area (not only irrigated areas)
    #~ inputFiles["irrigation_water_withdrawal"   ]    = inputDirectory + "/" + "irrigationWaterWithdrawal_annuaTot_output"
    #~ inputFiles["evaporation_from_irrigation"   ]    = inputDirectory + "/" + "evaporation_from_irrigation_annuaTot_output"
    #~ inputFiles["precipitation_at_irrigation"   ]    = inputDirectory + "/" + "precipitation_at_irrigation_annuaTot_output"

    # TODO: add info about area equipped with irrigation 
    # TODO: add info about return flow fraction for domestic, industry and livestock  

    # output that will be calculated 
    output = {}
    variable_names  = inputFiles.keys()
    #~ variable_names += ['irrigation_water_consumption']
    for var in variable_names:
        output[var] = {}
        output[var]['file_name'] = outputDirectory + "/" + str(var) + "_annual_country.nc"
        output[var]['unit']      = "km3.year-1"
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
    attributeDictionary['title'      ]  = "PCR-GLOBWB 2"
    attributeDictionary['institution']  = "Dept. of Physical Geography, Utrecht University"
    attributeDictionary['source'     ]  = "None"
    attributeDictionary['history'    ]  = "None"
    attributeDictionary['references' ]  = "None"
    attributeDictionary['comment'    ]  = "None"
    # additional attribute defined in PCR-GLOBWB 
    attributeDictionary['description'] = "prepared by Edwin H. Sutanudjaja"

    # initiate the netcd object: 
    tssNetCDF = MakingNetCDF(cloneMapFile = cloneMapFileName, \
                             attribute = attributeDictionary, \
                             cellSizeInArcMinutes = cellSizeInArcMinutes)
    # making netcdf files:
    for var in variable_names:
        tssNetCDF.createNetCDF(output[var]['file_name'], var, output[var]['unit'])

    # class (country) ids
    uniqueIDsFile = "/projects/0/dfguu/users/edwin/data/country_shp_from_tianyi/World_Polys_High.map"
    uniqueIDs = pcr.nominal(\
                vos.readPCRmapClone(uniqueIDsFile, cloneMapFileName, tmp_directory, 
                                    None, False, None, True))
    uniqueIDs = pcr.readmap(uniqueIDsFile)
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
        print fulldate

        # reading pcraster files:
        for var in inputFiles.keys():        
            
            # netcdf input file name:
            inputFile = inputFiles[var]
            if var!= "total_groundwater_abstraction":
                inputFile = inputFile + "_" + fulldate + "_to_" + fulldate + ".nc"
            if var in ["domesticGrossDemand",  \
                       "industryGrossDemand",  \
                       "livestockGrossDemand", \
                       "domesticNettoDemand",  \
                       "industryNettoDemand",  \
                       "livestockNettoDemand"]:
               inputFile = inputFiles[var]
            print inputFile   

            # fixing the year for water demand years (only between 1958 and 2010)
            fulldate_for_reading_netcdf = fulldate
            if iYear < 1960: fulldate_for_reading_netcdf = '%4i-%02i-%02i'  %(int(1960), int(12), int(31))
            if iYear > 2010: fulldate_for_reading_netcdf = '%4i-%02i-%02i'  %(int(2010), int(12), int(31))
            
            # reading PCR-GLOBWB values
            output[var]['pcr_value'] = vos.netcdf2PCRobjClone(ncFile = inputFile,\
                                                              varName = "Automatic",\
                                                              dateInput = fulldate_for_reading_netcdf,
                                                              useDoy = None,
                                                              cloneMapFileName  = cloneMapFileName,
                                                              LatitudeLongitude = True,
                                                              specificFillValue = None)
            # - water demand files/values are still in m.day-1 ; we have to convert them to m.year-1
            if var in ["domesticGrossDemand",  \
                       "industryGrossDemand",  \
                       "livestockGrossDemand", \
                       "domesticNettoDemand",  \
                       "industryNettoDemand",  \
                       "livestockNettoDemand"]:
                number_of_days_in_the_year = 365
                if calendar.isleap(iYear): number_of_days_in_the_year = 366 
                output[var]['pcr_value'] = output[var]['pcr_value'] * number_of_days_in_the_year
                                                              

        #~ # calculating irrigation water consumption
        #~ output['irrigation_water_consumption']['pcr_value'] = output['evaporation_from_irrigation']['pcr_value'] * \
                                                              #~ vos.getValDivZero(output['irrigation_water_withdrawal']['pcr_value'], \
                                                                                #~ output['irrigation_water_withdrawal']['pcr_value'] +\
                                                                                #~ output['precipitation_at_irrigation']['pcr_value'])
        
        # upscaling to the class (country) units and writing to netcdf files and a table
        for var in output.keys():
            
            print var
            
            # covering the map with zero
            pcrValue = pcr.cover(output[var]['pcr_value'], 0.0)
            
            # convert values from m to m3
            pcrValue =  pcrValue * cellArea

            # upscaling to the class (country) units and converting the units to km3/year
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
        txt_file = open(table_directory + "/" + "summary_" + fulldate + ".txt", "w")
        for var in output.keys():
            header += " " + str(var)
            header += "_km3"
            cmd    += " " + str(tmp_directory) + "/" + str(var) + ".tmp"
        cmd += " " + str(tmp_directory) + "/" + "summary_" + fulldate + ".txt.tmp"
        print cmd
        os.system(cmd)
        # - add header to txt file
        header += "\n" 
        txt_file.write(header)
        # - add map2col output to the txt_file
        map2col_file = open(tmp_directory+"/" + "summary_" + fulldate + ".txt.tmp", "r")
        txt_file.write(map2col_file.read())
        # - close all open txt files
        txt_file.close()
        map2col_file.close()
        
        # remove all temporary files
        cmd = 'rm -r '+ tmp_directory + "/*"
        os.system(cmd)
        
