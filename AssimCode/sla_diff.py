#!/public/home/lfs/miniforge3/envs/LFS-EnOI/bin/python

import xarray as xr
import sys
import numpy as np
def main():
    date_case = sys.argv[1]
    date = sys.argv[2]
    ttype = sys.argv[3]
    print(date),print(ttype)
    model_filename = "../modeldata/z0-"+str(date)[:4]+"-"+str(date)[4:6]+"-"+str(date)[6:]+"_25km.nc"
    model_filename1 = "../../Run_data_3days/Data_25km/output_"+str(date_case)+"_25km/z0-"+str(date)[:4]+"-"+str(date)[4:6]+"-"+str(date)[6:]+"_25km.nc"
    obs_filename   = "../input/SLA/"+str(date)[:4]+"/sla"+str(date)+".nc"
    ds_obs   = xr.open_dataset(obs_filename,decode_times=False)
    meanssh = xr.open_dataset('./ssh_ncdata/mdt_25km_m0.nc',decode_times=False)['mdt'].squeeze()
  #  ssh_trend = xr.open_dataset('./ssh_ncdata/sla_trend.nc',decode_times=False)
  #  ssh_trend0 = ssh_trend['sla_trend_93-08'].squeeze()
  #  ssh_trend1 = ssh_trend['sla_trend_09-24'].squeeze()
  #  ssh_season = xr.open_dataset('./ssh_ncdata/sla_season.nc')['sla_season'].sel(time='2001'+str(date)[4:]).squeeze()
    try:
        ds_model = xr.open_dataset(model_filename1,decode_times=False)
        print('try z_b in forecast')
        print(model_filename1)
    except Exception as e:
        try:
            print('try z_b in modeldata')
            ds_model = xr.open_dataset(model_filename,decode_times=False)
        except Exception as e:
            return print("Error opening dataset:", e)
    if ttype=='SLA':
        ds_model['ssh'] = ds_model['ssh'] - meanssh.data
 #   if (int(str(date)[:4])-2009) >=0:
 #       nyear0 = 5
 #       nyear1 = int(str(date)[:4])-2009
 #   else:
 #       nyear0 = int(str(date)[:4])-2004
 #       nyear1 = 0
 #   print('add seasonal change and trend')
 #   sla_season_trend = ssh_season.data + nyear0*ssh_trend0.data + nyear1*ssh_trend1.data
 #   ds_model['ssh'] = ds_model['ssh'] + sla_season_trend

    _, lat2d_model = np.meshgrid(ds_model['lon'].values, ds_model['lat'].values)
    area_model = np.cos(np.deg2rad(lat2d_model))
    ds_model = ds_model.assign(area=(('lat', 'lon'), area_model))    
    _, lat2d_obs = np.meshgrid(ds_obs['longitude'].values, ds_obs['latitude'].values)
    area_obs = np.cos(np.deg2rad(lat2d_obs))
    ds_obs = ds_obs.assign(area=(('latitude', 'longitude'), area_obs))
    avg_model = float(ds_model['ssh'].weighted(ds_model['area']).mean(dim=["lat", "lon"]))
    if ttype=='SLA':
        avg_obs   = float(ds_obs['sla'].weighted(ds_obs['area']).mean(dim=["latitude", "longitude"]))
    #    diff_2d = ds_model['ssh'].squeeze().data - ds_obs['sla'].squeeze().data
    else:
        avg_obs   = float(ds_obs['ssh'].weighted(ds_obs['area']).mean(dim=["latitude", "longitude"]))
     #   diff_2d = ds_model['ssh'].squeeze().data - ds_obs['ssh'].squeeze().data
 #   print('save results')

    diff = -(avg_obs - avg_model)
     
    print("SLA/SSH:model-obs=%.4f"%diff)
 #   diff_2d = -1*sla_season_trend+diff
 #   ds_out = xr.Dataset(
 #   {
 #       "sla_diff": (("lat", "lon"),diff_2d),
 #   },
 #   coords={
 #       "lat": ds_model['lat'].data,
 #       "lon": ds_model['lon'].data,
 #   }
 #   )
 #   ds_out.to_netcdf('sla_diff_'+str(date)+'.nc')
 #   diff_2d[np.isnan(diff_2d)] = 1e35
 #   flat = diff_2d.flatten()
 #   output_file = 'sla_diff.uf'
 #   with open(output_file, 'w') as f:
 #       for i in range(0, flat.size, 5):
 #           slice_ = flat[i:i+5]
 #           line = " ".join("{:E}".format(x) for x in slice_)
 #           f.write(line + "\n")
    with open("sla_diff", "w") as f:
        f.write(f"{diff:.6f}\n") 

if __name__ == '__main__':
    main()
