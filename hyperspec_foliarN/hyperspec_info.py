#add the necessary meta-data to the topo/BRDF corrected .h5 files 

#python hyperspec_info.py
import numpy as np
import h5py as h5
from os import listdir
from os.path import join


read_dir = 'Z:/SRER/Martha/hyperspec/flightlines/indv'
write_dir = 'Z:/SRER/Martha/hyperspec/flightlines/output'

for file in listdir(read_dir):
    if not file.endswith('.h5'):
        continue
    hyper_read = h5.File(join(read_dir,file), 'r')
    hyper_write = h5.File(join(write_dir,file), 'a')
    
    print(file)
    #get the site ID
    parts = file.split("_")
    siteID = parts[2]
    
    #create the group hierarcy in the output dataset
    hyper_write.create_group('SRER/Reflectance/Metadata/Spectral_Data')
    
    #get the wavelength dataset
    wavelength = hyper_read['SRER/Reflectance/Metadata/Spectral_Data/Wavelength']
    #copy the wavelength to the hdf write file
    hyper_write.create_dataset('SRER/Reflectance/Metadata/Spectral_Data/Wavelength', data=wavelength)
    
    #get the mapinfo
    mapinfo = hyper_read['SRER/Reflectance/Metadata/Coordinate_System/Map_Info']
    hyper_write.create_dataset('SRER/Reflectance/Metadata/Coordinate_System/Map_Info', data=mapinfo)
    
    #get attributes of the reflectance dataset
    refldata_read = hyper_read['SRER/Reflectance/Reflectance_Data']
    spatial_extent = refldata_read.attrs["Spatial_Extent_meters"]
    
    refldata_write = hyper_write['BRDF/Correction']
    refldata_write.attrs["Spatial_Extent_meters"] = spatial_extent
       
    #Extract bad band windows
    refl_read = hyper_read['SRER/Reflectance']
    badband1 = (refl_read.attrs["Band_Window_1_Nanometers"])
    badband2 = (refl_read.attrs["Band_Window_2_Nanometers"])
    
    refl_write = hyper_write['SRER/Reflectance']
    refl_write.attrs["Band_Window_1_Nanometers"] = badband1
    refl_write.attrs["Band_Window_2_Nanometers"] = badband2
    
    #Extract projection information
    proj = hyper_read['SRER/Reflectance/Metadata/Coordinate_System/Proj4']
    hyper_write.create_dataset('SRER/Reflectance/Metadata/Coordinate_System/Proj4', data=proj)
    
    epsg = hyper_read['SRER/Reflectance/Metadata/Coordinate_System/EPSG Code']
    hyper_write.create_dataset('SRER/Reflectance/Metadata/Coordinate_System/EPSG Code', data=epsg)
    
    hyper_write.close
    hyper_read.close