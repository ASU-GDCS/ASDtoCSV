# ASDtoCSV
Python code to convert an ASD spectrometer output file to a CSV format for easy uptake.  Output will be two columns with first column the wavelength and the second column the specified spectral output (default: Reflectance)

# Dependencies
Should work on any Python interpreter > Python 3.6 (only tested on 3.9.1).  This script uses only Python standard libraries.

# Usage
```
usage: asd2csv.py [-h] [--type {Reflectance,Raw,Reference}] [--sigdig SIGDIG] [--force_data_format {0,1,2,3}]
                  [--no_range_errors] [--dn]
                  ASD CSV

positional arguments:
  ASD                   Binary ASD file
  CSV                   Output file - will contain a column/entry for each valid asd file

optional arguments:
  -h, --help            show this help message and exit
  --type {Reflectance,Raw,Reference}, -t {Reflectance,Raw,Reference}
                        Which data to output, default=Reflectance
  --sigdig SIGDIG, -d SIGDIG
                        Format output to a specified number of significant digits
  --force_data_format {0,1,2,3}
                        Override header data format
  --no_range_errors     Ignore if some records have data outside specified dynamic range
  --dn                  Return DN values without normalization for 'Raw' data
```
