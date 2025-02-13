
    out_directory     = "/scratch-shared/edwinbar/electricity_water_demand_estimate_based_on_aqueduct_2021/historical/"
    # ~ out_directory = sys.argv[1]

    # output file code (which will be used as part of output file names)												
    output_file_code = "electricity_water_demand_historical"
    # ~ output_file_code = str(sys.argv[2])

    # industrial gross demand (m.day-1, monthly resolution)
    directory_for_water_demand   = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/historical_and_ssp_files/"
    industrial_gross_demand_file = directory_for_water_demand + "/" + "industry_water_demand_historical_1960-2019.nc"
    # ~ industrial_gross_demand  = directory_for_water_demand + "/" + str(sys.argv[3])

    # ~ staYear = 1960
    # ~ endYear = 2019
    staYear = int(sys.argv[4])
    endYear = int(sys.argv[5])

#~ edwinbar@fcn73.local.snellius.surf.nl:/home/edwinbar$ ls -lah /projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/historical_and_ssp_files/
#~ total 614G
#~ drwxr-xr-x 2 edwin edwin 4.0K Jan 21 16:19 .
#~ drwxr-xr-x 5 edwin edwin 4.0K Apr  1  2022 ..
#~ -rw-r--r-- 1 edwin edwin  51G Apr  1  2022 industry_water_demand_historical_1960-2019.nc
#~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 industry_water_demand_ssp1_2000-2100.nc
#~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 industry_water_demand_ssp3_2000-2100.nc
#~ -rw-r--r-- 1 edwin edwin  85G Apr  1  2022 industry_water_demand_ssp5_2000-2100.nc

out_directory = "/scratch-shared/edwinbar/electricity_water_demand_estimate_based_on_aqueduct_2021/historical/"
output_file_code = "electricity_water_demand_historical_1960-2019"
industrial_gross_demand_file = "industry_water_demand_historical_1960-2019.nc"
python zonal_calculation_tss_for_prometheus.py ${out_directory} ${output_file_code} ${industrial_gross_demand_file} 1960 2019 &

out_directory = "/scratch-shared/edwinbar/electricity_water_demand_estimate_based_on_aqueduct_2021/ssp1/"
output_file_code = "electricity_water_demand_ssp1_2000-2100"
industrial_gross_demand_file = "industry_water_demand_ssp1_2000-2100.nc"
python zonal_calculation_tss_for_prometheus.py ${out_directory} ${output_file_code} ${industrial_gross_demand_file} 2000 2100 &

out_directory = "/scratch-shared/edwinbar/electricity_water_demand_estimate_based_on_aqueduct_2021/ssp3/"
output_file_code = "electricity_water_demand_ssp3_2000-2100"
industrial_gross_demand_file = "industry_water_demand_ssp3_2000-2100.nc"
python zonal_calculation_tss_for_prometheus.py ${out_directory} ${output_file_code} ${industrial_gross_demand_file} 2000 2100 &

out_directory = "/scratch-shared/edwinbar/electricity_water_demand_estimate_based_on_aqueduct_2021/ssp5/"
output_file_code = "electricity_water_demand_ssp5_2000-2100"
industrial_gross_demand_file = "industry_water_demand_ssp5_2000-2100.nc"
python zonal_calculation_tss_for_prometheus.py ${out_directory} ${output_file_code} ${industrial_gross_demand_file} 2000 2100 &

wait

