#!/usr/bin/env python
import struct
import time
import argparse
import os
import glob
import math
from collections import OrderedDict

data_type_dict = {
    0:"Raw",
    1:"Reflectance",
    2:"Radiance",
    3:"No_Units",
    4:"Irradiance",
    5:"QI",
    6:"Transmittance",
    7:"Unknown",
    8:"Absorbance"
}

type_abbrev_dict = {
    "Raw":"Raw",
    "Reflectance":"Refl",
    "Reference":"Rfnc"
}

data_format_dict = {
    0:"numeric", 
    1:"integer", 
    2:"double", 
    3:"Unknown"
}

instrument_dict = {
    0:"Unknown", 
    1:"PSII",
    2:"LSVNIR",
    3:"FSVNIR",
    4:"FSFR",
    5:"FSNIR",
    6:"CHEM",
    7:"FSFR Unattended"
}

flag1_dict = {
    1:"nir saturation",
    2:"swir1 satruation",
    3:"swir2 saturation",
    8:"Tec1 alarm",
    16:"Tec2 alarm"
}

def dot(x,y):
    return sum([xi*yi for xi, yi in zip(x,y)])

def clean_envihdr_array(chunks, keyword, fmt=""):
    itemfmt = "{{:{}}}".format(fmt) if fmt != "" else "{}"
    str = "{} = {{ ".format(keyword)
    padstr = " "*len(str)
    str += ", ".join([itemfmt.format(itm) for itm in chunks[0]])
    for chunk in chunks[1:]:
        str += ",\n{}".format(padstr) + ", ".join([itemfmt.format(itm) for itm in chunk])
    str += "}\n"
    return str


class ASDSpec(object):
    """Read binary ASD file as a data frame"""
    def __init__(self, asdfile, range_errors=True, data_format=None,
                  do_normalize=True):
        self.filename = asdfile
        self.range_errors = range_errors
        ##Try to open file as a binary file
        bref = open(asdfile,"rb")
        ##Seek past signature
        bref.seek(3,0)

        ## Comments
        self.comments = bref.read(157).decode().replace('\x00',"")

        # Spectrum acquisition time
        bref.seek(182,0)
        seconds = struct.unpack('<L', bref.read(4))[0]
        self.aquisition_time = time.ctime(seconds)

        # Program and file version
        bref.seek(178,0)
        pv = struct.unpack("<B",bref.read(1))[0]
        self.program_version = "{}.{}".format(pv >> 4, pv & 7)
        fv = struct.unpack("<B",bref.read(1))[0] # idem for file version
        self.file_version = "{}.{}".format(fv >> 4, fv & 7)

        # Read the VNIR dark subtraction field.
        bref.seek(181,0)
        DC = struct.unpack("<B",bref.read(1))[0]
        if DC == 1:
            self.vnir_dark_sub = True
        else:
            if DC == 0:
                self.vnir_dark_sub = False
            else:
                self.vnir_dark_sub = NA
        # Read the dark spectrum datetime. The date and time are represented as the number of seconds since midnight on 1st
        # January 1970.
        seconds = struct.unpack("<L",bref.read(4))[0]
        self.dark_meas_time = time.ctime(seconds)

        # Read the spectrum data type. The type code is in range 0-8.

        bref.seek(186,0)
        self.datatype = data_type_dict[struct.unpack("<B",bref.read(1))[0]]

        # Read the reference spectrum datetime.
        seconds = struct.unpack("<L",bref.read(4))[0]
        self.white_meas_time = time.ctime(seconds)

        # Read GPS data.
        bref.seek(334,0)
        self.gps_trueHeading = struct.unpack("<d",bref.read(8))[0]
        self.gps_speed = struct.unpack("<d",bref.read(8))[0]
        self.gps_latitude = struct.unpack("<d",bref.read(8))[0]
        self.gps_longitude = struct.unpack("<d",bref.read(8))[0]
        self.gps_altitude = struct.unpack("<d",bref.read(8))[0]

        # Read the integration time.
        bref.seek(390,0)
        self.vnir_int_time = struct.unpack("<L",bref.read(4))[0]
        self.vnir_int_time_units = "ms"

        # Read the fore optic information.
        bref.seek(394,0)
        self.fore_optic = struct.unpack("<h",bref.read(2))[0]

        # Read the dark current correction value.
        bref.seek(396,0)
        self.dark_curr_corr = struct.unpack("<h",bref.read(2))[0]

        # Read the instrument number
        bref.seek(400,0)
        self.serial_number = str(struct.unpack("<h",bref.read(2))[0])

        # Read the warning flags
        bref.seek(421,0)
        warning_flags = struct.unpack("<4B",bref.read(4))
        #  if (sum(warningFlags)) {
        #    warning(paste("There appears to be a warning flag in the file:", f, "\nThis may indicate a problem with one of the detectors caused either by saturation (too much light) or by a failure of thermoelectric cooling."))
        #  }
        self.warning1 = "None"
        self.warning2 = "None"
        if warning_flags[0]:
            self.warning1 = "AVGFIXed"
        if warning_flags[1]:
            self.warning2 = flag1_dict[warning_flags[1]]

        # Read averaging information
        bref.seek(425,0)
        self.dark_curr_averaging = struct.unpack("<h",bref.read(2))[0] ##Num of DC measurements in the avg
        self.white_ref_averaging = struct.unpack("<h",bref.read(2))[0] ##Num of WR in the average
        self.averaging = struct.unpack("<h",bref.read(2))[0] ##Num of spec samples in the avg

        # Read the instrument model LS stands for LabSpec, FS for FieldSpec, FR for Full Range
        bref.seek(431,0)
        self.instrument_model = instrument_dict[struct.unpack("<B",bref.read(1))[0]]

        # Read the SWIR detector gain and offset settings.
        bref.seek(436,0)
        self.swir1_gain = struct.unpack("<h",bref.read(2))[0]
        self.swir2_gain = struct.unpack("<h",bref.read(2))[0]
        self.swir1_offset = struct.unpack("<h",bref.read(2))[0]
        self.swir2_offset = struct.unpack("<h",bref.read(2))[0]

        # Read the detector join wavelengths.
        bref.seek(444,0)
        self.join1_wavelength = struct.unpack("<f",bref.read(4))[0]
        self.join1_wavelength_units = "nm"
        self.join2_wavelength = struct.unpack("<f",bref.read(4))[0]
        self.join2_wavelength_units = "nm"

        # Read the smart detector data
        bref.seek(452,0)
        self.smart_detector_data = "/".join([str(sdd) for sdd in struct.unpack("<8f",bref.read(32))])
        #  if (sum(self.smart_detector_data))
        #    warning(paste("There appears to be data from a smart detector in the file:", f, "\nThis function does not support importing smart detector data"))
        #  }

        # Read spectra data. First reads some relevent information from the file header, then builds the wavelength scale and
        # reads the spectrum data values. If a reference spectrum is also present it will read that too.

        # Read the number of channels on the detector.
        bref.seek(204,0)
        self.num_channels = struct.unpack("<h",bref.read(2))[0]

        #  if self.num_channels == 2151:
        #      sensortype = "ASD FieldSpec Pro"
        #  elif self.num_channels == 751:
        #      sensortype = "ASD HH/VNIR"
        #  else:
        #      sensortype = "Unknown"

        # Read the wavelength information.
        bref.seek(191,0)
        wstart = struct.unpack("<f",bref.read(4))[0] # The first wavelengt
        bref.seek(195,0)
        wstep = struct.unpack("<f",bref.read(4))[0] # The interval between wavelengths.
        wend = wstart + self.num_channels * wstep - 1 # Calculate the last wavelengt

        # Build the wavelength scale
        self.wavelength = [step * wstep + wstart for step in range(self.num_channels)]
        self.wavelength_units = "nm"

        # Read the data format
        bref.seek(199,0)
        dataFormatCode = struct.unpack("<B",bref.read(1))[0] # In range 0 to 3.
        if data_format is not None:
            self.data_format = data_format_dict[data_format] # Format for arg in readBin
            struct_format = ["<{}f", "<{}l", "<{}d", None][data_format].format(self.num_channels)
            recbytes = [4, 4, 8, 0][data_format]
        else:
            self.data_format = data_format_dict[dataFormatCode] # Format for arg in readBin
            struct_format = ["<{}f", "<{}l", "<{}d", None][dataFormatCode].format(self.num_channels)
            recbytes = [4, 4, 8, 0][dataFormatCode]
        if recbytes < 1:
            raise RuntimeError("ASD records are of unknown datatype")

        # Read the instrument's dynamic range. This will be used for some basic validation.
        bref.seek(418,0)
        self.inst_dynamic_range = 2**struct.unpack("<h",bref.read(2))[0]

        # Read the target spectrum.  The 'Indico Version 8 File Format' document specifies that the spectrum starts at byte
        # 485. However it appears to actually start at byte 484.
        bref.seek(484,0)
        # The file format appears to have changed with file version and even file pre-processing (raw and ref) The following
        # code guess the size argument based on the number of channels it should retrieve
        self.rawdata = struct.unpack(struct_format,bref.read(recbytes*self.num_channels))

        ##Get the reference data if provided
        self.filesize = bref.seek(0,2)
        self.refdesc = ""
        if (484+2*recbytes*self.num_channels) < self.filesize:
            bref.seek(484+recbytes*self.num_channels)
            refbool, refsec, specsec = struct.unpack("<Hqq",bref.read(18))
            refdescsize = struct.unpack("<h",bref.read(2))[0]
            self.refdesc = bref.read(refdescsize).decode().replace('\x00',"")

            #bref.seek(-recbytes*self.num_channels,2)
            self.refdata = list(struct.unpack(struct_format,bref.read(recbytes*self.num_channels)))
        else:
            self.refdata = None

        # If any target spectrum data values lie outside the dynamic range of the instrument this probably indicates that
        # something has gone wrong. This could be due to an incorrect offset or data type when reading the binary file.
        if self.range_errors and sum([d>self.inst_dynamic_range for d in self.rawdata]) > 0:
            raise RuntimeError("ASD records are larger than dynamic range: "+\
                    "{}".format(2**self.inst_dynamic_range))

        # Normalize the target spectrum
        def normalize(wl, val):
            if wl <= self.join1_wavelength:
                return val / self.vnir_int_time
            if (wl > self.join1_wavelength) and (wl <= self.join2_wavelength):
                return val * self.swir1_gain / 2048
            if wl > self.join2_wavelength:
                return val * self.swir2_gain / 2048
            return val

        if do_normalize and (self.datatype == "Raw"):
            self.rawdata = [normalize(self.wavelength[i], self.rawdata[i]) for i in range(self.num_channels)]

    def transform(self, output_type="Reflectance"):
        if output_type == "Reflectance":
            if self.refdata is not None:
                return [self.rawdata[i] / self.refdata[i] for i in range(self.num_channels)]
            elif self.datatype == "Reflectance":
                return list(self.rawdata)
            else:
                raise RuntimeError("Reference data not contained in file")
        elif output_type == "Raw":
            return list(self.rawdata)
        elif output_type == "Reference":
            if self.refdata is None:
                raise RuntimeError("Reference data not contained in file")
            return list(self.refdata)
        elif output_type == "Transmittance":
            if self.datatype == "Transmittance":
                return list(self.rawdata)
            else:
                raise RuntimeError(f"Cannot create Transmittance spectra from header type {self.datatype}")
        else:
            raise RuntimeError("Output type {} not understood".format(output_type))

    def to_csv(self, csvname, output_type="Reflectance", value_fmt=""):
        vals = self.transform(output_type)
        ##Name the wavelength field,
        wlname = "wl_{}".format(self.wavelength_units)
        ##Name spectral field
        abbrev = type_abbrev_dict[output_type]
        if abbrev != "":
            colname = os.path.splitext(os.path.basename(self.filename))[0]+"_"+abbrev
        else:
            colname = os.path.splitext(os.path.basename(self.filename))[0]
        ##Start writing
        with open(csvname, 'w') as oref:
            #Write header
            oref.write(f"{wlname},{colname}\n")
            for row in range(len(self.wavelength)):
                oref.write(f"{self.wavelength[row]},{value_fmt.format(vals[row])}\n")
        return

def main(args):
    if not os.path.exists(args.input):
        raise RuntimeError(f"Path '{args.input}' not found")

    try:
        asdf = ASDSpec(args.input,range_errors=not args.no_range_errors,data_format=args.data_format,do_normalize=args.dn)
    except Exception as exc:
        raise RuntimeError("ERROR - cannot input data from file {} - {}".format(args.input,str(exc)))

    ##Write odict to requested format
    print("Writing data to file {}".format(args.output))
    valfmt = "{{{}}}".format(":.{}f".format(args.sigdig) if args.sigdig else "")
    asdf.to_csv(args.output, output_type=args.type, value_fmt=valfmt)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--type", "-t", default="Reflectance", choices=["Reflectance","Raw","Reference"], help="Which data to output, default=Reflectance")
    parser.add_argument("--sigdig", "-d", default=None, type=int, help="Format output to a specified number of significant digits")
    parser.add_argument("--force_data_format", dest="data_format", default=None, type=int, choices=list(data_format_dict.keys()), help="Override header data format")
    parser.add_argument("--no_range_errors", action="store_true", help="Ignore if some records have data outside specified dynamic range")
    parser.add_argument("--dn", action="store_false", help="Return DN values without normalization for 'Raw' data")
    parser.add_argument("input", metavar="ASDFILE/DIR", help="Binary ASD file or directory containing asd files")
    parser.add_argument("output", metavar="FILE", help="Output file - will contain a column/entry for each valid asd file")
    args = parser.parse_args()
    main(args)
