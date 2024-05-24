rm(list=ls(all=T))
library('SPEI')
library(ncdf4) # package for netcdf manipulation
#install.packages('raster')
library(raster) # package for raster manipulation
#install.packages('rgdal')
#install.packages("maps")
library(maps)
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
setwd('~/iCloud/com~apple~CloudDocs/uconn/classes/enve5810_hydrometeorology/project/raw_data/')

################################################################################

# declaring initial var retrieval functions and very inefficiently pulling data
# netcdf files

lon_getter <- function(path_to_clt) {
  ins_clt_out = nc_open(path_to_clt)
  ins_lon = ncvar_get(ins_clt_out, "lon")
  nc_close(ins_clt_out)
  return(ins_lon)
}

lat_getter <- function(path_to_clt) {
  ins_clt_out = nc_open(path_to_clt)
  ins_lat = ncvar_get(ins_clt_out, "lat")
  nc_close(ins_clt_out)
  return(ins_lat)
}

var_getter <- function(path_from_raw, var_name) {
  var_out = nc_open(path_from_raw)
  var_array = ncvar_get(var_out, var_name)
  fillValue = ncatt_get(var_out, var_name, "_FillValue")
  var_array[var_array==fillValue$value] = NA
  nc_close(var_out)
  var_vec_long = as.matrix(var_array)
  return(var_vec_long)
}


bcc_clt_vec1 = var_getter("./output2/clt_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc", "clt")
bcc_huss_vec1 = var_getter("./output2/huss_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc", "huss")
bcc_pr_vec1 = var_getter("./output2/pr_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc", "pr")
bcc_ps_vec1 = var_getter("./output2/ps_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc", "ps")
bcc_sfcWind_vec1 = var_getter("./output2/sfcWind_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc", "sfcWind")
bcc_tasmax_vec1 = var_getter("./output2/tasmax_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc", "tasmax")
bcc_tasmin_vec1 = var_getter("./output2/tasmin_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc", "tasmin")
bcc_clt_vec1


bcc_clt_vec2 = var_getter("./output2/clt_Amon_BCC-CSM2-MR_207101-210012_uscrop_2.nc", "clt")
bcc_huss_vec2 = var_getter("./output2/huss_Amon_BCC-CSM2-MR_207101-210012_uscrop_2.nc", "huss")
bcc_pr_vec2 = var_getter("./output2/pr_Amon_BCC-CSM2-MR_207101-210012_uscrop_2.nc", "pr")
bcc_ps_vec2 = var_getter("./output2/ps_Amon_BCC-CSM2-MR_207101-210012_uscrop_2.nc", "ps")
bcc_sfcWind_vec2 = var_getter("./output2/sfcWind_Amon_BCC-CSM2-MR_207101-210012_uscrop_2.nc", "sfcWind")
bcc_tasmax_vec2 = var_getter("./output2/tasmax_Amon_BCC-CSM2-MR_207101-210012_uscrop_2.nc", "tasmax")
bcc_tasmin_vec2 = var_getter("./output2/tasmin_Amon_BCC-CSM2-MR_207101-210012_uscrop_2.nc", "tasmin")

can_clt_vec1 = var_getter("./output2/clt_Amon_CanESM5_201501-204412_uscrop_2.nc", "clt")
can_huss_vec1 = var_getter("./output2/huss_Amon_CanESM5_201501-204412_uscrop_2.nc", "huss")
can_pr_vec1 = var_getter("./output2/pr_Amon_CanESM5_201501-204412_uscrop_2.nc", "pr")
can_ps_vec1 = var_getter("./output2/ps_Amon_CanESM5_201501-204412_uscrop_2.nc", "ps")
can_sfcWind_vec1 = var_getter("./output2/sfcWind_Amon_CanESM5_201501-204412_uscrop_2.nc", "sfcWind")
can_tasmax_vec1 = var_getter("./output2/tasmax_Amon_CanESM5_201501-204412_uscrop_2.nc", "tasmax")
can_tasmin_vec1 = var_getter("./output2/tasmin_Amon_CanESM5_201501-204412_uscrop_2.nc", "tasmin")

can_clt_vec2 = var_getter("./output2/clt_Amon_CanESM5_207101-210012_uscrop_2.nc", "clt")
can_huss_vec2 = var_getter("./output2/huss_Amon_CanESM5_207101-210012_uscrop_2.nc", "huss")
can_pr_vec2 = var_getter("./output2/pr_Amon_CanESM5_207101-210012_uscrop_2.nc", "pr")
can_ps_vec2 = var_getter("./output2/ps_Amon_CanESM5_207101-210012_uscrop_2.nc", "ps")
can_sfcWind_vec2 = var_getter("./output2/sfcWind_Amon_CanESM5_207101-210012_uscrop_2.nc", "sfcWind")
can_tasmax_vec2 = var_getter("./output2/tasmax_Amon_CanESM5_207101-210012_uscrop_2.nc", "tasmax")
can_tasmin_vec2 = var_getter("./output2/tasmin_Amon_CanESM5_207101-210012_uscrop_2.nc", "tasmin")

miroc_clt_vec1 = var_getter("./output2/clt_Amon_MIROC6_201501-204412_uscrop_2.nc", "clt")
miroc_huss_vec1 = var_getter("./output2/huss_Amon_MIROC6_201501-204412_uscrop_2.nc", "huss")
miroc_pr_vec1 = var_getter("./output2/pr_Amon_MIROC6_201501-204412_uscrop_2.nc", "pr")
miroc_ps_vec1 = var_getter("./output2/ps_Amon_MIROC6_201501-204412_uscrop_2.nc", "ps")
miroc_sfcWind_vec1 = var_getter("./output2/sfcWind_Amon_MIROC6_201501-204412_uscrop_2.nc", "sfcWind")
miroc_tasmax_vec1 = var_getter("./output2/tasmax_Amon_MIROC6_201501-204412_uscrop_2.nc", "tasmax")
miroc_tasmin_vec1 = var_getter("./output2/tasmin_Amon_MIROC6_201501-204412_uscrop_2.nc", "tasmin")

miroc_clt_vec2 = var_getter("./output2/clt_Amon_MIROC6_207101-210012_uscrop_2.nc", "clt")
miroc_huss_vec2 = var_getter("./output2/huss_Amon_MIROC6_207101-210012_uscrop_2.nc", "huss")
miroc_pr_vec2 = var_getter("./output2/pr_Amon_MIROC6_207101-210012_uscrop_2.nc", "pr")
miroc_ps_vec2 = var_getter("./output2/ps_Amon_MIROC6_207101-210012_uscrop_2.nc", "ps")
miroc_sfcWind_vec2 = var_getter("./output2/sfcWind_Amon_MIROC6_207101-210012_uscrop_2.nc", "sfcWind")
miroc_tasmax_vec2 = var_getter("./output2/tasmax_Amon_MIROC6_207101-210012_uscrop_2.nc", "tasmax")
miroc_tasmin_vec2 = var_getter("./output2/tasmin_Amon_MIROC6_207101-210012_uscrop_2.nc", "tasmin")

mri_clt_vec1 = var_getter("./output2/clt_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc", "clt")
mri_huss_vec1 = var_getter("./output2/huss_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc", "huss")
mri_pr_vec1 = var_getter("./output2/pr_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc", "pr")
mri_ps_vec1 = var_getter("./output2/ps_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc", "ps")
mri_sfcWind_vec1 = var_getter("./output2/sfcWind_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc", "sfcWind")
mri_tasmax_vec1 = var_getter("./output2/tasmax_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc", "tasmax")
mri_tasmin_vec1 = var_getter("./output2/tasmin_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc", "tasmin")

mri_clt_vec2 = var_getter("./output2/clt_Amon_MRI-ESM2-0_207101-210012_uscrop_2.nc", "clt")
mri_huss_vec2 = var_getter("./output2/huss_Amon_MRI-ESM2-0_207101-210012_uscrop_2.nc", "huss")
mri_pr_vec2 = var_getter("./output2/pr_Amon_MRI-ESM2-0_207101-210012_uscrop_2.nc", "pr")
mri_ps_vec2 = var_getter("./output2/ps_Amon_MRI-ESM2-0_207101-210012_uscrop_2.nc", "ps")
mri_sfcWind_vec2 = var_getter("./output2/sfcWind_Amon_MRI-ESM2-0_207101-210012_uscrop_2.nc", "sfcWind")
mri_tasmax_vec2 = var_getter("./output2/tasmax_Amon_MRI-ESM2-0_207101-210012_uscrop_2.nc", "tasmax")
mri_tasmin_vec2 = var_getter("./output2/tasmin_Amon_MRI-ESM2-0_207101-210012_uscrop_2.nc", "tasmin")

################################################################################

# doing some basic conversions into units required for Penman-Monteith

bcc_sfcWind_2m_vec1 = bcc_sfcWind_vec1 * 0.7479511
bcc_sfcWind_2m_vec2 = bcc_sfcWind_vec2 * 0.7479511

can_sfcWind_2m_vec1 = can_sfcWind_vec1 * 0.7479511
can_sfcWind_2m_vec2 = can_sfcWind_vec2 * 0.7479511

miroc_sfcWind_2m_vec1 = miroc_sfcWind_vec1 * 0.7479511
miroc_sfcWind_2m_vec2 = miroc_sfcWind_vec2 * 0.7479511

mri_sfcWind_2m_vec1 = mri_sfcWind_vec1 * 0.7479511
mri_sfcWind_2m_vec2 = mri_sfcWind_vec2 * 0.7479511

month_list = c(29.5, 29.5, 30.5, 30.5, 30.5, 30.5, 31, 30.5, 30.5, 30.5, 30.5, 31)

bcc_pr_mm_vec1 = bcc_pr_vec1 * 86400 * month_list
bcc_pr_mm_vec2 = bcc_pr_vec2 * 86400 * month_list

can_pr_mm_vec1 = can_pr_vec1 * 86400 * month_list
can_pr_mm_vec2 = can_pr_vec2 * 86400 * month_list

miroc_pr_mm_vec1 = miroc_pr_vec1 * 86400 * month_list
miroc_pr_mm_vec2 = miroc_pr_vec2 * 86400 * month_list

mri_pr_mm_vec1 = mri_pr_vec1 * 86400 * month_list
mri_pr_mm_vec2 = mri_pr_vec2 * 86400 * month_list

bcc_tasmin_degC_vec1 = bcc_tasmin_vec1 - 273.15
bcc_tasmin_degC_vec2 = bcc_tasmin_vec2 - 273.15

can_tasmin_degC_vec1 = can_tasmin_vec1 - 273.15
can_tasmin_degC_vec2 = can_tasmin_vec2 - 273.15

miroc_tasmin_degC_vec1 = miroc_tasmin_vec1 - 273.15
miroc_tasmin_degC_vec2 = miroc_tasmin_vec2 - 273.15

mri_tasmin_degC_vec1 = mri_tasmin_vec1 - 273.15
mri_tasmin_degC_vec2 = mri_tasmin_vec2 - 273.15

bcc_tasmax_degC_vec1 = bcc_tasmax_vec1 - 273.15
bcc_tasmax_degC_vec2 = bcc_tasmax_vec2 - 273.15

can_tasmax_degC_vec1 = can_tasmax_vec1 - 273.15
can_tasmax_degC_vec2 = can_tasmax_vec2 - 273.15

miroc_tasmax_degC_vec1 = miroc_tasmax_vec1 - 273.15
miroc_tasmax_degC_vec2 = miroc_tasmax_vec2 - 273.15

mri_tasmax_degC_vec1 = mri_tasmax_vec1 - 273.15
mri_tasmax_degC_vec2 = mri_tasmax_vec2 - 273.15

bcc_ps_Pa_vec1 = bcc_ps_vec1 * 0.001
bcc_ps_Pa_vec2 = bcc_ps_vec2 * 0.001

can_ps_Pa_vec1 = can_ps_vec1 * 0.001
can_ps_Pa_vec2 = can_ps_vec2 * 0.001

miroc_ps_Pa_vec1 = miroc_ps_vec1 * 0.001
miroc_ps_Pa_vec2 = miroc_ps_vec2 * 0.001

mri_ps_Pa_vec1 = mri_ps_vec1 * 0.001
mri_ps_Pa_vec2 = mri_ps_vec2 * 0.001

bcc_lon = lon_getter("./output2/clt_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc")
bcc_lat = lat_getter("./output2/clt_Amon_BCC-CSM2-MR_201501-204412_uscrop_2.nc")

can_lon = lon_getter("./output2/clt_Amon_CanESM5_201501-204412_uscrop_2.nc")
can_lat = lat_getter("./output2/clt_Amon_CanESM5_201501-204412_uscrop_2.nc")

miroc_lon = lon_getter("./output2/clt_Amon_MIROC6_201501-204412_uscrop_2.nc")
miroc_lat = lat_getter("./output2/clt_Amon_MIROC6_201501-204412_uscrop_2.nc")

mri_lon = lon_getter("./output2/clt_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc")
mri_lat = lat_getter("./output2/clt_Amon_MRI-ESM2-0_201501-204412_uscrop_2.nc")

time_obs1 = seq(as.Date("2015-01-16"), by = "month", length.out = 360)
time_obs2 = seq(as.Date("2071-01-16"), by = "month", length.out = 360)

mon_seq1 = 1:360
mon_seq2 = 841:1200

#bcc_coord1 = as.matrix(expand.grid(bcc_lon, bcc_lat, mon_seq1))
#bcc_coord2 = as.matrix(expand.grid(bcc_lon, bcc_lat, mon_seq2))

dim(bcc_clt_vec1)
bcc_clt_vec1
combo = expand.grid(bcc_lon, bcc_lat)
bcc_clt_arr1 = array(bcc_clt_vec1, dim = c(85, 49, 360))
bcc_clt_arr1[,,2]


# combines individual variable vectors into an array of all input variables
arr_comb = function(clt_vec, pr_vec, ps_vec, sfcWind_vec, tasmin_vec, tasmax_vec, nlon, nlat) {
  arr_out = array(dim = c(nlon, nlat, 360, 6))
  var_list = list(clt_vec, pr_vec, ps_vec, sfcWind_vec, tasmin_vec, tasmax_vec)
  i = 1
  for (var_vec in var_list) {
    arr_temp = array(var_vec, dim = c(nlon, nlat, 360))
    arr_out[,,,i] = arr_temp
    i = i + 1
  }
  return(arr_out)
}

bcc_2044_arr = arr_comb(bcc_clt_vec1, bcc_pr_mm_vec1, bcc_ps_Pa_vec1, bcc_sfcWind_2m_vec1, bcc_tasmin_degC_vec1, bcc_tasmax_degC_vec1, 85, 49)
bcc_2100_arr = arr_comb(bcc_clt_vec2, bcc_pr_mm_vec2, bcc_ps_Pa_vec2, bcc_sfcWind_2m_vec2, bcc_tasmin_degC_vec2, bcc_tasmax_degC_vec2, 85, 49)

can_2044_arr = arr_comb(can_clt_vec1, can_pr_mm_vec1, can_ps_Pa_vec1, can_sfcWind_2m_vec1, can_tasmin_degC_vec1, can_tasmax_degC_vec1, 34, 19)
can_2100_arr = arr_comb(can_clt_vec2, can_pr_mm_vec2, can_ps_Pa_vec2, can_sfcWind_2m_vec2, can_tasmin_degC_vec2, can_tasmax_degC_vec2, 34, 19)

miroc_2044_arr = arr_comb(miroc_clt_vec1, miroc_pr_mm_vec1, miroc_ps_Pa_vec1, miroc_sfcWind_2m_vec1, miroc_tasmin_degC_vec1, miroc_tasmax_degC_vec1, 68, 39)
miroc_2100_arr = arr_comb(miroc_clt_vec2, miroc_pr_mm_vec2, miroc_ps_Pa_vec2, miroc_sfcWind_2m_vec2, miroc_tasmin_degC_vec2, miroc_tasmax_degC_vec2, 68, 39)

mri_2044_arr = arr_comb(mri_clt_vec1, mri_pr_mm_vec1, mri_ps_Pa_vec1, mri_sfcWind_2m_vec1, mri_tasmin_degC_vec1, mri_tasmax_degC_vec1, 85, 49)
mri_2100_arr = arr_comb(mri_clt_vec2, mri_pr_mm_vec2, mri_ps_Pa_vec2, mri_sfcWind_2m_vec2, mri_tasmin_degC_vec2, mri_tasmax_degC_vec2, 85, 49)


# calculates Penman-Monteith PET and outputs a new array with PET and CWB added on
pen_calc2 = function(work_arr, work_lat) {
  nlon = dim(work_arr)[1]
  nlat = dim(work_arr)[2]
  arr_out = array(dim = c(nlon, nlat, 360, 8))
  for (lat_int in 1:nlat) {
    for (lon_int in 1:nlon) {
      #arr_temp = array(dim = c(nlon, nlat, 360, 8))
      #print(lat_int)
      #print(lon_int)
      pen = array(penman(Tmin = work_arr[lon_int, lat_int, , 5], Tmax = work_arr[lon_int, lat_int, , 6], U2 = work_arr[lon_int, lat_int, , 4], lat = work_lat[lat_int], CC = work_arr[lon_int, lat_int, , 1], P = work_arr[lon_int, lat_int, , 3], verbose = F))
      arr_out[lon_int, lat_int, , 7] = pen
    }
  }
  #arr_out[,,,8] = arr_out[,,,2] - arr_out[,,,7]
  arr_out[,,,c(-7, -8)] = work_arr
  arr_out[,,,8] = arr_out[,,,2] - arr_out[,,,7]
  return(arr_out)
}


# calculates 1 and 12-month SPEI and outputs a new array with these values added on 
spei_calc2 = function(work_arr, end_year) {
  nlon = dim(work_arr)[1]
  nlat = dim(work_arr)[2]
  arr_out = array(dim = c(nlon, nlat, 360, 10))
  arr_out[,,,c(-9,-10)] = work_arr
  for (lat_int in 1:nlat) {
    print(lat_int/nlat*100)
    for (lon_int in 1:nlon) {
      ts_temp = ts(work_arr[lon_int, lat_int, , 8], end = c(end_year, 12), frequency = 12)
      spei1 = spei(ts_temp, 1, verbose = F)
      spei12 = spei(ts_temp, 12, verbose = F)
      for (time_int in 1:360) {
        arr_out[lon_int, lat_int, time_int, 9] = spei1$fitted[time_int]
        arr_out[lon_int, lat_int, time_int, 10] = spei12$fitted[time_int]
      }
    }
  }
  return(arr_out)
}

bcc_2044_pen_arr = pen_calc2(bcc_2044_arr, bcc_lat)
bcc_2100_pen_arr = pen_calc2(bcc_2100_arr, bcc_lat)

can_2044_pen_arr = pen_calc2(can_2044_arr, can_lat)
can_2100_pen_arr = pen_calc2(can_2100_arr, can_lat)

miroc_2044_pen_arr = pen_calc2(miroc_2044_arr, miroc_lat)
miroc_2100_pen_arr = pen_calc2(miroc_2100_arr, miroc_lat)

mri_2044_pen_arr = pen_calc2(mri_2044_arr, mri_lat)
mri_2100_pen_arr = pen_calc2(mri_2100_arr, mri_lat)



bcc_2044_full_arr = spei_calc2(bcc_2044_pen_arr, 2044)
bcc_2100_full_arr = spei_calc2(bcc_2100_pen_arr, 2100)

can_2044_full_arr = spei_calc2(can_2044_pen_arr, 2044)
can_2100_full_arr = spei_calc2(can_2100_pen_arr, 2100)

miroc_2044_full_arr = spei_calc2(miroc_2044_pen_arr, 2044)
miroc_2100_full_arr = spei_calc2(miroc_2100_pen_arr, 2100)

mri_2044_full_arr = spei_calc2(mri_2044_pen_arr, 2044)
mri_2100_full_arr = spei_calc2(mri_2100_pen_arr, 2100)

# writes netcdf file, combining 1 and 12 month scale calculations into single file
write_netcdf = function(work_arr, end_year, ins_lon, ins_lat, path_including_file) {
  work_arr[is.infinite(work_arr)] = NA
  if (end_year == 2044) time_seq = 1:360
  if (end_year == 2100) time_seq = 841:1200
  nlon = dim(work_arr)[1]
  nlat = dim(work_arr)[2]
  londim = ncdim_def("lon", units = "degrees_east", ins_lon)
  latdim = ncdim_def("lat", units = "degrees_north", ins_lat)
  timedim = ncdim_def("time", units = "month", longname = "number of months since December 16, 2014", vals = time_seq)
  
  fillvalue = 1e32
  
  clt_def = ncvar_def("clt", units = "%", list(londim, latdim, timedim), fillvalue, longname = "average cloud coverage", prec = "double")
  pr_def = ncvar_def("pr", units = "mm", list(londim, latdim, timedim), fillvalue, longname = "average precipitation", prec = "double")
  ps_def = ncvar_def("ps", units = "kPa", list(londim, latdim, timedim), fillvalue, longname = "average surface pressure", prec = "double")
  sfcWind_def = ncvar_def("sfcWind", units = "m/s", list(londim, latdim, timedim), fillvalue, longname = "average wind speed at 2 m above surface", prec = "double")
  tasmin_def = ncvar_def("tasmin", units = "celsius", list(londim, latdim, timedim), fillvalue, longname = "average minimum daily surface temperature", prec = "double")
  tasmax_def = ncvar_def("tasmax", units = "celsius", list(londim, latdim, timedim), fillvalue, longname = "average maximum daily surface temperature", prec = "double")
  pet_def = ncvar_def("pet", units = "mm", list(londim, latdim, timedim), fillvalue, longname = "Penman-Monteith calculated monthly potential evapotranspiration", prec = "double")
  bal_def = ncvar_def("bal", units = "mm", list(londim, latdim, timedim), fillvalue, longname = "monthly water balance", prec = "double")
  speivar1_def = ncvar_def("spei1mon", units = "", list(londim, latdim, timedim), fillvalue, longname = "1 month SPEI", prec = "double")
  speivar12_def = ncvar_def("spei12mon", units = "", list(londim, latdim, timedim), fillvalue, longname = "12 month SPEI", prec = "double")
  
  ncout = nc_create(path_including_file, list(clt_def, pr_def, ps_def, sfcWind_def, tasmin_def, tasmax_def, pet_def, bal_def, speivar1_def, speivar12_def), force_v4 = T)
  ncvar_put(ncout, clt_def, work_arr[,,,1])
  ncvar_put(ncout, pr_def, work_arr[,,,2])
  ncvar_put(ncout, ps_def, work_arr[,,,3])
  ncvar_put(ncout, sfcWind_def, work_arr[,,,4])
  ncvar_put(ncout, tasmin_def, work_arr[,,,5])
  ncvar_put(ncout, tasmax_def, work_arr[,,,6])
  ncvar_put(ncout, pet_def, work_arr[,,,7])
  ncvar_put(ncout, bal_def, work_arr[,,,8])
  ncvar_put(ncout, speivar1_def, work_arr[,,,9])
  ncvar_put(ncout, speivar12_def, work_arr[,,,10])
  
  ncatt_put(ncout, "lon", "axis", "X")
  ncatt_put(ncout, "lat", "axis", "Y")
  ncatt_put(ncout, "time", "axis", "T")
}


# write arrays to netcdf files, must uncomment lines to complete file write
# write_netcdf(bcc_2044_full_arr, 2044, bcc_lon, bcc_lat, "./spei_data/bcc_2044_spei.nc")
# write_netcdf(bcc_2100_full_arr, 2100, bcc_lon, bcc_lat, "./spei_data/bcc_2100_spei.nc")
# 
# write_netcdf(can_2044_full_arr, 2044, can_lon, can_lat, "./spei_data/can_2044_spei.nc")
# write_netcdf(can_2100_full_arr, 2100, can_lon, can_lat, "./spei_data/can_2100_spei.nc")
# 
# write_netcdf(miroc_2044_full_arr, 2044, miroc_lon, miroc_lat, "./spei_data/miroc_2044_spei.nc")
# write_netcdf(miroc_2100_full_arr, 2100, miroc_lon, miroc_lat, "./spei_data/miroc_2100_spei.nc")
# 
# write_netcdf(mri_2044_full_arr, 2044, mri_lon, mri_lat, "./spei_data/mri_2044_spei.nc")
# write_netcdf(mri_2100_full_arr, 2100, mri_lon, mri_lat, "./spei_data/mri_2100_spei.nc")
