from setuptools import setup, Extension
import setuptools
import shutil, os.path, subprocess

def call_pkgconfig(args):
    ans = subprocess.check_output(args)
    return ans.split()
cflags = ['-std=c++14','-DPICOJSON_USE_INT64']
cflags.append('-I/home/james/work/oqtcpp/include')
cflags.append('-fvisibility=hidden')
cflags.append('-Wshadow')
libs =['-L/usr/local/lib', '-loqt','-lgeos_c']

srcs = ['src/oqttiles.cpp', 'src/geos_wrapper.cpp', 'src/prepare_geometries.cpp', 'src/mvt.cpp', 'src/maketiledata.cpp']

from distutils.command.build_ext import build_ext
from distutils.sysconfig import customize_compiler

class my_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        build_ext.build_extensions(self)

ext_modules = [
    Extension(
        'oqttiles._oqttiles',
        srcs,
        include_dirs=[
            #'/home/james/build/boost_1_65_1',
            '/usr/local/include',
            '/home/james/work/oqtcpp/include',
            #'/home/james/build/mapbox_geometry/include',
            
        ],
        extra_link_args=libs,
        extra_compile_args=cflags
        
    ),
    
]



setup(
    name='oqttiles',
    packages=['oqttiles'],
    version='0.0.1',
    long_description='',
    ext_modules=ext_modules,
    zip_safe=False,
    cmdclass = {'build_ext': my_build_ext}
)
