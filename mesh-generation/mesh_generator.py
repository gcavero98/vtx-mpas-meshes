import argparse
import os

import xarray as xr

from personalized_variable_resolution import variable_resolution_latlonmap, \
    view_resolution_map

from jigsaw_generator import jigsaw_gen_sph_grid

from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.mesh.conversion import convert
from mpas_tools.io import write_netcdf


def get_mesh_from_resolution(resolution_ds, basename='./mesh'):
    print('\n>> Generating an MPAS mesh')

    # jigsaw
    print('\n\t .- Jigsaw generation')
    mesh_file = jigsaw_gen_sph_grid(resolution_ds['resolution'].values,
                                    resolution_ds['lon'].values,
                                    resolution_ds['lat'].values,
                                    basename=basename)

    # mpas-tools

    print('\n\t .- Jigsaw to netCDF')
    out_file_triangles = basename + '.triangles.nc'
    jigsaw_to_netcdf(msh_filename=mesh_file,
                     output_name=out_file_triangles,
                     on_sphere=True, sphere_radius=1.0)

    print('\n\t .- Convert to MPAS format')
    out_file_mpas = basename + '.grid.nc'
    write_netcdf(
        convert(xr.open_dataset(out_file_triangles),
                dir=os.path.dirname(basename),
                graphInfoFileName=basename + ".graph.info"),
        out_file_mpas)

    return out_file_mpas


def full_generation_process(mpas_grid_file, grid, redo=True, **kwargs):
    if os.path.isfile(mpas_grid_file) and not redo:
        print(' .. already available')
        return
    os.system('rm -f ' + mpas_grid_file)

    resolution_ds = variable_resolution_latlonmap(grid, **kwargs)

    os.system('rm -rf tmp/')
    os.system('mkdir -p tmp/')

    tmp_mesh_file = get_mesh_from_resolution(resolution_ds,
                                             basename='tmp/mesh')

    os.system('cp ' + tmp_mesh_file + ' ' + mpas_grid_file)
    os.system('rm -rf tmp/')
    if not os.path.isfile(mpas_grid_file):
        raise IOError('The file we had to generate was not generated')

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "-g", "--grid", type=str, default='doughnut',
        help="""
            Grid option: \n 
            "  doughnut": 
                     <High resolution> area of a certain radius <size>. 
                     Linear increase of resolution to a certain <low
                     resolution> value after <margin>km. 
                     Keep the constant low resolution value for a while
                     (10*low_resolution)km and then increase it linearly
                     again until reaching 1000km (to save space).
                     The requested MPAS region should be circular and 
                     have a radius of <size>+<margin>. The buffer 
                     generated by the MPAS-Limited-Area code will then
                     consist of a few "rings" of <lowresolution> cells.
                     \n
            """
    )

    parser.add_argument(
        "-highr", "--highresolution", required=True, type=float,
        help="Highest-resolution of grid (km).",
    )

    parser.add_argument(
        "-lowr", "--lowresolution", default=None, type=float,
        help="Lowest-resolution of grid (km).",
    )

    parser.add_argument(
        "-size", "--size", default=None, type=float,
        help="Radius of the highest-resolution area of the grid (km).",
    )

    parser.add_argument(
        "-margin", "--margin", default=None, type=float,
        help="Size of the variable resolution boundary around the "
             "high resolution area (km).",
    )

    parser.add_argument(
        "-clat", "--clat", default=0., type=float,
        help="Central latitude.",
    )

    parser.add_argument(
        "-clon", "--clon", default=0., type=float,
        help="Central longitude.",
    )

    parser.add_argument(
        "-n", "--name", default="doughnut", type=str,
        help="output basename for directory and files."
    )

    # -p generates plots
    parser.add_argument(
        "-p", "--withplots", action="store_true",
        help="generate plots to view the mesh.",
    )

    # -o overwrite
    parser.add_argument(
        "-o", "--overwrite", action="store_true",
        help="overwrite existing folder.",
    )

    args = parser.parse_args()

    if args.name == '':
        raise ValueError('Please give a non trivial name.')

    folder = 'data/' + args.name + '/'
    if os.path.isdir(folder):
        if not args.overwrite:
            print('Sure to overwrite?')
            raise IOError('For security, overwriting is disabled. Give '
                          'different tests different names or erase the'
                          'existing folder: ' + folder)
        else:
            print('Overwriting folder ' + folder)
            os.system('rm -rf ' + folder)

    os.system('mkdir -p ' + folder)
    basename = folder + args.name

    ds = variable_resolution_latlonmap(args.grid,
                                       highresolution=args.highresolution,
                                       lowresolution=args.lowresolution,
                                       size=args.size,
                                       margin=args.margin,
                                       lat_ref=args.clat,
                                       lon_ref=args.clon,
                                       )

    print(ds)

    if args.withplots:
        print('Plotting')
        view_resolution_map(ds, pdfname=basename + '.resolution.pdf',
                            list_distances=[
                                            #1500,
                                            500,
                                            ds.attrs['border'],
                                            ds.attrs['radius']
                                            ])

    mesh_file = get_mesh_from_resolution(ds, basename=basename)

    print('DONE. This is the mesh ' + mesh_file)
    mpas_ds = xr.open_dataset(mesh_file)
    print(mpas_ds)
