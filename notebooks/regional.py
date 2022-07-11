#!/usr/bin/env python
import numpy as np
import os
import xarray as xr

# In[2]:

from vtxmpasmeshes.mesh_generator import variable_resolution_latlonmap

# define grid type and the desired mesh info
grid = 'doughnut'
test = {'highresolution': 3, 'lowresolution': 20, 'size': 30,
        'margin': 80, 'lat_ref': 40., 'lon_ref': 0.}

# we generate an Xarray Dataset containing, for each lat, lon point,
# the distance to a reference point and the corresponding resolution
# given the mesh info
if not os.path.exists('map.nc'):
    latlon_map = variable_resolution_latlonmap(grid, **test)
    latlon_map.to_netcdf('map.nc')
else:
    latlon_map = xr.open_dataset('map.nc')

# In[3]:


from vtxmpasmeshes.mesh_generator import get_mesh_from_resolution

# temporary folder
tmp_dir = 'tmp'
if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)

# we generate an spherical grid using a simplification of a Jigsaw script
# and we save it as a NetCDF file in the temporary folder
global_grid = f'{tmp_dir}/mesh.grid.nc'
if not os.path.exists(global_grid):
    global_grid, graph_info = get_mesh_from_resolution(
        resolution_ds=latlon_map, basename=f'{tmp_dir}/mesh'
    )

grid = xr.open_dataset(global_grid)

# In[4]:


from vtxmpasmeshes.mesh_generator import cut_circular_region

regional_grid = 'circle.grid.nc'
if not os.path.exists(f'grids/{regional_grid}'):
    grids_dir = 'grids'
    if not os.path.exists(grids_dir):
        os.makedirs(grids_dir)
    # we cut the global mesh given the radius (in meters) and the center point
    regional_grid, grid_info = cut_circular_region(
        mpas_global_file=global_grid,
        region_radius_meters=latlon_map.attrs['radius'] * 1000,
        lat_cen=latlon_map.attrs['lat_ref'],
        lon_cen=latlon_map.attrs['lon_ref']
    )
    # os.system('rm -rf tmp points.txt')
    os.system('mv circle* global* grids')

# In[5]:


from vtxmpasmeshes.mpas_plots import view_mpas_regional_mesh

# we can plot the regional MPAS mesh with a reanalysis grid.
# do_plot_resolution_rings and do_plot_wrf_grid must be set to False
view_mpas_regional_mesh(mpas_grid_file=f'grids/{regional_grid}',
                        do_plot_resolution_rings=False,
                        do_plot_era5_grid=True,
                        do_plot_wrf_grid=False,
                        **{'vmax': 3})
exit()
# In[6]:


######################################################################

# Open dataset and update attributes
regional_ds = xr.open_dataset(f'grids/{regional_grid}')
num_boundary_layers = 8

for name, value in latlon_map.attrs.items():
    regional_ds.attrs['vtx-param-' + str(name)] = value
print(regional_ds)
regional_ds.attrs['vtx-region-num_boundary_layers'] = num_boundary_layers
lowres = latlon_map.attrs['lowresolution']
radius = latlon_map.attrs['radius']
region_border = radius + (num_boundary_layers * lowres) * 0.6
regional_ds.attrs['vtx-region_border'] = region_border

regional_grid_2 = 'circle_2.grid.nc'
regional_ds.to_netcdf(f'grids/{regional_grid_2}')

print(regional_ds)


# In[9]:


latlon_map = view_mpas_regional_mesh(mpas_grid_file=f'grids/{regional_grid_2}',
                                     do_plot_resolution_rings=True,
                                     do_plot_era5_grid=True,
                                     do_plot_wrf_grid=True,
                                     **{'vmax': 7.})





