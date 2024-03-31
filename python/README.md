## Installation

A list of dependencies can be found in `requirements/base.in`. It is recommended to create a python environment:

```
python -m venv env
```

To load/unload the environment:

```
source env/bin/activate
deactivate
```

Using pip-compile (`pip install pip-tools`), compile a `.txt` file and install all dependencies:

```
cd requirements
pip-compile base.in                  #This can take a while
python -m pip-install -r base.txt
```

## Usage
