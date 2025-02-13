
import os
import sys
import shutil
import datetime

import numpy as np
import netCDF4 as nc
import pcraster as pcr

import virtualOS as vos

# General steps:
# ~ 0. irr_area_5min = min(cell_area_5min, irr_area_5min)
# ~ 1. downscaled_irr_area_30sec = (irr_area_30sec / sum_irr_area_30sec) * irr_area_5min
# ~ 2. downscaled_irr_area_30sec = min(cell_area_30sec, downscaled_irr_area_30sec)
# ~ 3. remaining_irr_area_5min   = max(0.0, irr_area_5min - sum_downscaled_irr_area_30sec)
# ~ 4. downscaled_irr_area_30sec = downscaled_irr_area_30sec + (not_irr_assigned_yet_cell_area_30sec / sum_not_irr_assigned_yet_cell_area_30sec) * remaining_irr_area_5min 


class MakingNetCDF():
    
    def __init__(self, cloneMapFile, attribute = None):
        		
        # cloneMap
        # - the cloneMap must be at 5 arc min resolution
        cloneMap = pcr.readmap(cloneMapFile)
        cloneMap = pcr.boolean(1.0)
        
        # latitudes and longitudes
        self.latitudes  = np.unique(pcr.pcr2numpy(pcr.ycoordinate(cloneMap), vos.MV))[::-1]
        self.longitudes = np.unique(pcr.pcr2numpy(pcr.xcoordinate(cloneMap), vos.MV))

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

        rootgrp= nc.Dataset(ncFileName, 'w', format = self.format)

        #-create dimensions - time is unlimited, others are fixed
        rootgrp.createDimension('time', None)
        rootgrp.createDimension('lat', len(self.latitudes))
        rootgrp.createDimension('lon', len(self.longitudes))

        date_time= rootgrp.createVariable('time', 'f4', ('time',))
        date_time.standard_name = 'time'
        date_time.long_name = 'Days since 1901-01-01'

        date_time.units = 'Days since 1901-01-01' 
        date_time.calendar = 'standard'

        lat= rootgrp.createVariable('lat','f4',('lat',))
        lat.long_name = 'latitude'
        lat.units = 'degrees_north'
        lat.standard_name = 'latitude'

        lon= rootgrp.createVariable('lon','f4',('lon',))
        lon.standard_name = 'longitude'
        lon.long_name = 'longitude'
        lon.units = 'degrees_east'

        lat[:]= self.latitudes
        lon[:]= self.longitudes

        shortVarName = varName
        var= rootgrp.createVariable(shortVarName,'f4', ('time','lat','lon',), fill_value = vos.MV, zlib = True)
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

        shortVarName = varName        

        date_time = rootgrp.variables['time']
        date_time[posCnt] = nc.date2num(timeStamp, date_time.units, date_time.calendar)

        rootgrp.variables[shortVarName][posCnt,:,:] = (varField)

        rootgrp.sync()
        rootgrp.close()

def main():


    # output and temporary directories
    out_directory     = "/scratch-shared/edwinbar/electricity_water_demand/test/"
    # ~ out_directory = sys.argv[1]
    tmp_directory     = out_directory + "/" + "tmp" + "/"
    # - making output and temporary directories
    if os.path.exists(out_directory):
        shutil.rmtree(out_directory)
    os.makedirs(out_directory)
    os.makedirs(tmp_directory)
    # - moving to the output directory
    os.chdir(out_directory)

    # - output table directory
    table_directory = out_directory + "/table/"
    os.makedirs(table_directory)


    # output file code (which will be used as part of output file names)												
    output_file_code = "electricity_water_demand_historical"
    # ~ output_file_code = str(sys.argv[2])


    # industrial gross demand (m.day-1, monthly resolution)
    directory_for_water_demand   = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/historical_and_ssp_files/"
    industrial_gross_demand_file = directory_for_water_demand + "/" + "industry_water_demand_historical_1960-2019.nc"
    # ~ industrial_gross_demand  = directory_for_water_demand + "/" + str(sys.argv[3])
    

    # clone map
    cloneMapFileName = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/lddsound_05min_version_20210330.map"
    pcr.setclone(cloneMapFileName) 
    

    # class (country) ids
    uniqueIDsFile = "/home/edwinbar/github/edwinkost/pgb_data_for_prometheus_primes/iso3_countries/row_number_CNTR_RG_01M_2020_4326.shp.tif.map"
    # class (country) ids
    uniqueIDs = pcr.nominal(\
                vos.readPCRmapClone(uniqueIDsFile, cloneMapFileName, tmp_directory, 
                                    None, False, None, True))
    uniqueIDs = pcr.ifthen(pcr.scalar(uniqueIDs) >= 0.0, uniqueIDs)
    

    # landmask                               
    landmask05minFile = cloneMapFileName
    landmask = pcr.defined(pcr.readmap(landmask05minFile))
    landmask = pcr.ifthen(landmask, landmask)
    # - extending landmask with uniqueIDs
    landmask = pcr.cover(landmask, pcr.defined(uniqueIDs))
    

    # cell area at 5 arcmin resolution (unit: m2)
    cell_area_5min_file  = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/cdo_gridarea_clone_global_05min_correct_lats.nc"
    cell_area_5min       = pcr.cover(vos.netcdf2PCRobjCloneWithoutTime(ncFile  = cell_area_5min_file, \
                                                                       varName = "automatic", cloneMapFileName = cloneMapFileName, LatitudeLongitude = True, specificFillValue = None, absolutePath = None), 0.0)
    cellArea = pcr.ifthen(landmask, cell_area_5min)


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


    # get a sample cell for every id
    x_min_for_each_id = pcr.areaminimum(pcr.xcoordinate(pcr.boolean(1.0)), uniqueIDs)
    sample_cells      = pcr.xcoordinate(pcr.boolean(1.0)) == x_min_for_each_id
    y_min_for_each_id = pcr.areaminimum(pcr.ycoordinate(sample_cells), uniqueIDs)
    sample_cells      = pcr.ycoordinate(sample_cells) == y_min_for_each_id
    uniqueIDs_sample  = pcr.ifthen(sample_cells, uniqueIDs)
    # - save it to a pcraster map file
    pcr.report(uniqueIDs_sample, "sample.ids")                                


    # table contains country values of fraction of electricy water demand from industrial water demand
    # - column 1 is the row_number_CNTR_RG_01M_2020_4326.shp.tif.map
    # - column 2 is the fraction of electricy water demand from industrial water demand
    table_pow_man_split       = "/home/edwinbar/github/edwinkost/pgb_data_for_prometheus_primes/pow_man_split/pow_man_split.txt"
    # - country map of this fraction
    country_pow_man_split_map = pcr.lookupscalar(table_pow_man_split, uniqueIDs)
    # - cover missing values with 0.5
    country_pow_man_split_map = pcr.cover(country_pow_man_split_map, 0.5)
    # - make sure that every country has a unique value
    country_pow_man_split_map = pcr.areaaverage(country_pow_man_split_map, uniqueIDs)
    

    # start year and end year
    staYear = 2000
    endYear = 2005
    # ~ staYear = int(sys.argv[4])
    # ~ endYear = int(sys.argv[5])
    

    # attribute for netCDF files 
    attributeDictionary = {}
    attributeDictionary['title'      ]  = "Estimate of electricty water demand at the country scale."
    attributeDictionary['institution']  = "Dept. of Physical Geography, Utrecht University"
    attributeDictionary['source'     ]  = directory_for_water_demand
    attributeDictionary['history'    ]  = "None"
    attributeDictionary['references' ]  = "None"
    attributeDictionary['comment'    ]  = "The country code is based on the " +  str(uniqueIDsFile) + "."
    # additional attribute defined in PCR-GLOBWB 
    attributeDictionary['description']  = "Created by Edwin H. Sutanudjaja, see https://github.com/edwinkost/pgb_data_for_prometheus_primes/blob/main/scripts/zonal_calculation_tss_for_aqueduct_2021_for_electricity_water_demand.py"


    # initiate the netcd object: 
    tssNetCDF = MakingNetCDF(cloneMapFile = cloneMapFileName, \
                             attribute    = attributeDictionary)
    # initiate the netcdf file:
    output = {}
    var = "electricity_water_demand"
    output[var] = {}
    output[var]['file_name'] = "/scratch-shared/edwinbar/electricity_water_demand/test/" + output_file_code + ".nc"
    output[var]['unit']      = "m3.year-1"
    variable_names = [var]
    for var in variable_names:
        tssNetCDF.createNetCDF(output[var]['file_name'], var, output[var]['unit'])

    
    # calculating and writing to the netcdf file
    for iYear in range(staYear, endYear+1):
        
        print(iYear)
        
        for iMonth in range(0, 12+1):

            # initating the variable to calculate yearly electricity water demand
            if iMonth == 0: country_annual_electricity_water_demand_volume = 0.0

            # time stamp for reading netcdf file:
            fulldate = '%4i-%02i-%02i'  %(int(iYear), int(iMonth), int(1))

            # industrial water demand 
            # - unit: m.day-1
            industrial_water_demand = vos.netcdf2PCRobjClone(ncFile  = industrial_gross_demand_file,\
                                                             varName = "industryGrossDemand", dateInput = fulldate, useDoy = None, cloneMapFileName = cloneMapFileName, LatitudeLongitude = True, specificFillValue = None)
            # - number of days within a month
            number_of_days_in_this_month = calendar.monthrange(iYear, iMonth)[1]
            # - unit: m3.month-1
            monthly_industrial_water_demand_volume = industrial_water_demand * cellArea * number_of_days_in_this_month
            
            # country scale industrial water demand (unit: m3.month-1 per country)
            country_monthly_industrial_water_demand_volume = pcr.areatotal(monthly_industrial_water_demand_volume, uniqueIDs)
            
            # country scale electricity water demand (unit: m3.month-1 per country)
            country_monthly_electricity_water_demand_volume = country_pow_man_split_map * country_monthly_industrial_water_demand_volume
            
            # accumulate it over the year 
            country_annual_electricity_water_demand_volume  = country_annual_electricity_water_demand_volume + country_monthly_electricity_water_demand_volume


        # index for writing the netcdf
        if iYear == staYear: index = 0
        index = index + 1

        # write values to a netcdf file
        var = "electricity_water_demand"
        ncFileName = output[var]['file_name']
        varField = pcr.pcr2numpy(pcrValue, vos.MV)
        tssNetCDF.writePCR2NetCDF(ncFileName, var, varField, timeStamp, posCnt = index - 1)
            
        # plot the values at sample cells only and write values to a temporary pcraster map
        pcrFileName = str(tmp_directory) + "/" + str(var) + ".tmp"
        pcr.report(pcr.ifthen(pcr.defined(uniqueIDs_sample), country_annual_electricity_water_demand_volume), pcrFileName)

        
        # write class values to a table
        # - command line to call map2col
        cmd    = 'map2col -x 1 -y 2 -m NA sample.ids'
        # - header for the table
        header = "x y class_id"
        # - txt file that contains the table
        txt_file = open(table_directory + "/" + output_file_code + "_summary_" + str(iYear) + ".txt", "w")
        for var in output.keys():
            header += " " + str(var)
            cmd    += " " + str(tmp_directory) + "/" + str(var) + ".tmp"
        cmd += " " + str(tmp_directory) + "/" + output_file_code + "_summary_" + str(iYear) + ".txt.tmp"
        print(cmd)
        os.system(cmd)
        # - add header to txt file
        header += "\n" 
        txt_file.write(header)
        # - add map2col output to the txt_file
        map2col_file = open(tmp_directory + "/" + output_file_code + "_summary_" + str(iYear) + ".txt.tmp", "r")
        txt_file.write(map2col_file.read())
        # - close all open txt files
        txt_file.close()
        map2col_file.close()
        
        # remove all temporary files
        cmd = 'rm -r '+ tmp_directory + "/*"
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    sys.exit(main())


