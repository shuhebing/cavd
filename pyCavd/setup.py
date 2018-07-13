#!/usr/bin/env python3

#from distutils.core import setup
from setuptools import setup
from setuptools.extension import Extension
try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True


includedirs=["../voro++/src"]
libdirs = ["../voro++/src"]

#print use_cython
#use_cython = False
cmd_class = {}

if use_cython:
    ext = '.pyx'
    cmd_class.update({'build_ext':build_ext})
else:
    ext = '.cpp'

netstorage_srcfiles = [
        'zeo/netstorage'+ext, '../zeo++/networkstorage.cc', 
        '../zeo++/mindist.cc', '../zeo++/geometry.cc', '../zeo++/networkinfo.cc',
        '../zeo++/networkio.cc', '../zeo++/grid.cc', '../zeo++/symbcalc.cc',
        '../zeo++/string_additions.cc', 
        '../zeo++/voronoicell.cc', 
        '../zeo++/networkanalysis.cc', '../zeo++/graphstorage.cc', '../zeo++/area_and_volume.cc',
        '../zeo++/network.cc', '../zeo++/OMS.cc', '../zeo++/v_network.cc', '../zeo++/symmetry.cc',
        '../zeo++/networkaccessibility.cc', '../zeo++/channel.cc', '../zeo++/net.cc', '../zeo++/ray.cc',
	'../zeo++/rmsd.cc','../zeo++/material.cc', '../zeo++/psd.cc',
        ]
netinfo_srcfiles = ['zeo/netinfo'+ext, '../zeo++/networkinfo.cc']
netio_srcfiles = [
        'zeo/netio'+ext, '../zeo++/networkio.cc', 'zeo/netinfo'+ext, 
        #'../zeo++/networkinfo.cc', 'zeo/string_add.pxd', '../zeo++/string_additions.cc', 
        '../zeo++/networkinfo.cc',  '../zeo++/string_additions.cc', 
        '../zeo++/grid.cc', '../zeo++/mindist.cc', '../zeo++/symbcalc.cc',  '../zeo++/symmetry.cc',
        '../zeo++/networkstorage.cc', '../zeo++/geometry.cc', '../zeo++/net.cc', '../zeo++/rmsd.cc',
        ]
graphstorage_srcfiles = ['zeo/graphstorage'+ext, '../zeo++/graphstorage.cc']
psd_srcfiles = ['zeo/psd'+ext, '../zeo++/psd.cc']
voronoicell_srcfiles = [
        'zeo/voronoicell'+ext, '../zeo++/voronoicell.cc', '../zeo++/geometry.cc',
	'../zeo++/networkstorage.cc', '../zeo++/net.cc', '../zeo++/mindist.cc', 
	'../zeo++/networkinfo.cc', '../zeo++/rmsd.cc', '../zeo++/symmetry.cc', 
	'../zeo++/string_additions.cc', '../zeo++/ray.cc', '../zeo++/channel.cc', 
	'../zeo++/network.cc', '../zeo++/OMS.cc', '../zeo++/area_and_volume.cc', '../zeo++/networkaccessibility.cc', 
	'../zeo++/graphstorage.cc', '../zeo++/networkanalysis.cc', '../zeo++/v_network.cc',
        ]
channel_srcfiles = ['zeo/channel'+ext, '../zeo++/channel.cc']
highaccuracy_srcfiles = [
        'zeo/high_accuracy'+ext, '../zeo++/sphere_approx.cc', '../zeo++/networkstorage.cc', 
        '../zeo++/networkinfo.cc', '../zeo++/mindist.cc', '../zeo++/geometry.cc', '../zeo++/net.cc', 
	'../zeo++/symmetry.cc', '../zeo++/string_additions.cc', '../zeo++/ray.cc',
	'../zeo++/networkaccessibility.cc', '../zeo++/network.cc', '../zeo++/networkio.cc',
	'../zeo++/grid.cc', '../zeo++/symbcalc.cc', '../zeo++/voronoicell.cc', '../zeo++/graphstorage.cc',
	'../zeo++/channel.cc', '../zeo++/v_network.cc', '../zeo++/networkanalysis.cc',
	'../zeo++/area_and_volume.cc', '../zeo++/rmsd.cc', '../zeo++/material.cc', '../zeo++/psd.cc',
        ]
areavol_srcfiles = [
        'zeo/area_volume'+ext, '../zeo++/area_and_volume.cc', '../zeo++/networkinfo.cc', 
        '../zeo++/networkstorage.cc', '../zeo++/mindist.cc', '../zeo++/geometry.cc', 
        '../zeo++/networkio.cc', '../zeo++/grid.cc', '../zeo++/symbcalc.cc',
        '../zeo++/string_additions.cc', 'zeo/voronoicell'+ext, '../zeo++/voronoicell.cc', 
        '../zeo++/networkanalysis.cc', '../zeo++/graphstorage.cc', '../zeo++/symmetry.cc', 
        '../zeo++/network.cc', '../zeo++/OMS.cc', '../zeo++/v_network.cc', '../zeo++/ray.cc', '../zeo++/rmsd.cc',
        '../zeo++/networkaccessibility.cc', '../zeo++/channel.cc', '../zeo++/net.cc'
        ]
cluster_srcfiles = [
        'zeo/cluster'+ext, '../zeo++/cluster.cc', '../zeo++/networkstorage.cc',
        '../zeo++/networkinfo.cc', '../zeo++/mindist.cc', '../zeo++/geometry.cc', 
        '../zeo++/network.cc', '../zeo++/OMS.cc', '../zeo++/voronoicell.cc', '../zeo++/graphstorage.cc',
        '../zeo++/networkanalysis.cc', '../zeo++/channel.cc', '../zeo++/v_network.cc',
        '../zeo++/area_and_volume.cc',  '../zeo++/networkaccessibility.cc', 
        '../zeo++/string_additions.cc', '../zeo++/sphere_approx.cc', '../zeo++/net.cc',
	'../zeo++/symmetry.cc', '../zeo++/ray.cc', '../zeo++/rmsd.cc', '../zeo++/material.cc', '../zeo++/psd.cc',
        ]
cycle_srcfiles = [
        'zeo/cycle'+ext, '../zeo++/cycle.cc', '../zeo++/networkstorage.cc',
        '../zeo++/networkinfo.cc', '../zeo++/mindist.cc', '../zeo++/geometry.cc', 
        '../zeo++/network.cc', '../zeo++/OMS.cc', '../zeo++/voronoicell.cc', '../zeo++/graphstorage.cc',
        '../zeo++/networkanalysis.cc', '../zeo++/channel.cc', '../zeo++/v_network.cc',
        '../zeo++/area_and_volume.cc',  '../zeo++/networkaccessibility.cc', 
        '../zeo++/string_additions.cc', '../zeo++/sphere_approx.cc'
        ]
geometry_srcfiles = ['zeo/geometry'+ext, '../zeo++/geometry.cc']
setup(
    name = 'zeo',
    version = '0.1',
    description = "Python interface to Zeo++",
    url = "http://www.maciejharanczyk.info/Zeopp/",
    author = "Bharat Medasani, Maciej Haranzcyk",
    author_email = "bkmedasani@lbl.gov",
    packages=["zeo",],
    package_dir={"zeo": "zeo"},
    license = "",
    #cmdclass = {'build_ext':build_ext},
    cmdclass = cmd_class,
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Cython",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering"
        ],
    ext_modules = [
                   Extension("zeo.netstorage",
                             sources=netstorage_srcfiles, 
                             include_dirs=includedirs,
                             libraries = ['voro++'], 
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.geometry", 
                             sources=geometry_srcfiles,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.netinfo", 
                             sources=netinfo_srcfiles,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.voronoicell", 
                             sources=voronoicell_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'], 
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.netio",
                             sources=netio_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.graphstorage",
                             sources=graphstorage_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.psd", 
                             sources=psd_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.channel", 
                             sources=channel_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.high_accuracy", 
                             sources=highaccuracy_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.area_volume", 
                             sources=areavol_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.cluster", 
                             sources=cluster_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                   Extension("zeo.cycle", 
                             sources=cycle_srcfiles,
                             include_dirs=includedirs,
                             libraries = ['voro++'],
                             library_dirs = libdirs,
                             extra_compiler_args = ['-fPIC'],
                             language = 'c++'),
                             ]
)
