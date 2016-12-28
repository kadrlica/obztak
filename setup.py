from setuptools import setup, find_packages
import versioneer
import glob

NAME = 'obztak'
CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
Programming Language :: Python
Natural Language :: English
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Astronomy
Operating System :: MacOS
Operating System :: POSIX
License :: OSI Approved :: MIT License
"""
URL = 'https://github.com/obztak/%s'%NAME
DESC = "Bizarro observation tactician for DECam"
LONG_DESC = "See %s"%URL

setup(
    name=NAME,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url=URL,
    author='Alex Drlica-Wagner',
    author_email='kadrlica@fnal.gov',
    scripts = glob.glob('bin/*'),
    install_requires=[
        'numpy >= 1.7',
        'scipy >= 0.10.1',
        'basemap >= 1.0.6',
        'pandas >= 0.8.1',
        'setuptools',
    ],
    packages=find_packages(),
    package_data={'maglites':['data/*.dat','data/*.txt','data/*.csv']},
    description=DESC,
    long_description=LONG_DESC,
    platforms='any',
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f]
)
