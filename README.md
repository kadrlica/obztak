# Bizarro ObzTak

[![Build](https://img.shields.io/travis/kadrlica/obztak.svg)](https://travis-ci.org/kadrlica/obztak)
[![PyPI](https://img.shields.io/pypi/v/obztak.svg)](https://pypi.python.org/pypi/obztak)
[![Release](https://img.shields.io/github/release/kadrlica/obztak.svg)](../../releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](../../)

Bizarro observation tactician for the Dark Energy Camera (DECam) at CTIO. Used for the **Mag**ellanic Satel**lite**s **S**urvey (MagLiteS) and the **BL**anco **I**maging of the **S**outhern **S**ky (BLISS), and the **DE**Cam **L**ocal **V**olume Survey (DELVE).

### Installation

Bizarro ObzTak is now distributed via `pip`:
```
pip install obztak
```

### Running

The first step is to run `prepare_survey` to set up the necessary survey characterization files. Specifically, this script builds a list of survey fields and a list of expected time windows. These files will be written to the current directory by default.
```
> ./bin/prepare_survey --help
usage: prepare_survey [-h] [-p] [-f FIELDS] [-w WINDOWS]

Decide which fields to observe and time windows to observe.

optional arguments:
  -h, --help            show this help message and exit
  -p, --plot            Plot output. (default: False)
  -f FIELDS, --fields FIELDS
                        List of all target fields. (default:
                        target_fields.txt)
  -w WINDOWS, --windows WINDOWS
                        List of observation windows. (default:
                        observation_windows.txt)
```


There are two pimary executables: `survey_simulator` and `survey_observer`. Both use the same underlying architecture, but the `survey_simulator` simulates the entire survey while `survey_observer` is used to simulate only a specific chuck of the survey and create an output json file. Both have the `-p` option for real-time plotting.
```
> ./bin/survey_simulator --help
usage: survey_simulator [-h] [-v] [-p] [-fields FIELDS] [-w WINDOWS] [-d DONE]
                        [-o OUTFILE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Output verbosity.
  -p, --plot            Plot output.
  -fields FIELDS, --fields FIELDS
                        List of all target fields.
  -w WINDOWS, --windows WINDOWS
                        List of observation windows.
  -d DONE, --done DONE  List of fields that have been observed.
  -o OUTFILE, --outfile OUTFILE
                        Save output fields surveyed.
```

```
usage: survey_observer [-h] [-v] [-p] [-fields FIELDS] [-w WINDOWS] [-d DONE]
                       [-o OUTFILE] [--tstart TSTART] [--tstop TSTOP]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Output verbosity.
  -p, --plot            Plot output.
  -fields FIELDS, --fields FIELDS
                        List of all target fields.
  -w WINDOWS, --windows WINDOWS
                        List of observation windows.
  -d DONE, --done DONE  List of fields that have been observed.
  -o OUTFILE, --outfile OUTFILE
                        Save output fields surveyed.
  --tstart TSTART       Start time for observation.
  --tstop TSTOP         Stop time for observation.
```

## Development

To clone repository go to the directory where you want to work
```
git clone https://github.com/kadrlica/obztak.git
```
To keep your local copy up to date, pull changes from the remote repository to your local copy
```
git pull --all 
```
To commit an update
```
# git add command moves changes from the working directory to the staging area
git add example.txt 
# git commit takes the staged snapshot and commits it to the project history
git commit -m 'my comments' 
# git push moves a local branch or series of commits to main repository
git push -u origin master 
```
To work on a specific branch
```
# Example with the maglites branch...
git fetch origin maglites
# setup a new branch to track origin/maglites
git checkout -b maglites origin/maglites 
# push maglites branch (which is tracking origin/maglites) back to GitHub
git push origin maglites 
```
For more details on branches, see [here](http://stackoverflow.com/q/1783405/4075339).
