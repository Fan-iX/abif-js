/**
* Parse abi data
* 
* ABIF files are a proprietary binary sanger sequencing chromatogram data file
* created by Applied Biosystems (see
* https://projects.nfstc.org/workshops/resources/articles/ABIF_File_Format.pdf
* ). Data is parsed into an object containing all fields in the ABIF file.
* Fields may vary by machine and basecaller versions.
* This function is based on the `read.abif` function from R package `sangerseqR`.
* @param {ArrayBuffer} rawdata
* @returns {{header, directory, data}}
*/
function parse_abif(rawdata) {
    var enc = new TextDecoder("utf-8");
    function RTC(ab) {
        return enc.decode(ab);
    }
    function UInt8(ab) {
        var n = ab.byteLength
        if (n == 1) {
            return new DataView(ab).getUint8();
        }
        var arr = new Array(n);
        var dview = new DataView(ab);
        for (var i = 0; i < n; i++) {
            arr[i] = dview.getUint8(i);
        }
        return arr
    }

    function UInt16(ab) {
        var n = ab.byteLength / 2
        if (n == 1) {
            return new DataView(ab).getUint16();
        }
        var arr = new Array(n);
        var dview = new DataView(ab);
        for (var i = 0; i < n; i++) {
            arr[i] = dview.getUint16(i * 2);
        }
        return arr
    }

    function SInt16(ab) {
        var n = ab.byteLength / 2
        if (n == 1) {
            return new DataView(ab).getInt16();
        }
        var arr = new Array(n);
        var dview = new DataView(ab);
        for (var i = 0; i < n; i++) {
            arr[i] = dview.getInt16(i * 2);
        }
        return arr
    }

    function SInt32(ab) {
        var n = ab.byteLength / 4
        if (n == 1) {
            return new DataView(ab).getInt32();
        }
        var arr = new Array(n);
        var dview = new DataView(ab);
        for (var i = 0; i < n; i++) {
            arr[i] = dview.getInt32(i * 4);
        }
        return arr
    }

    function f32(ab) {
        var n = ab.byteLength / 4
        if (n == 1) {
            return new DataView(ab).getFloat32(0, true);
        }
        var arr = new Array(n);
        var dview = new DataView(ab);
        for (var i = 0; i < n; i++) {
            arr[i] = dview.getFloat32(i * 4, true);
        }
        return arr
    }

    function f64(ab) {
        var n = ab.byteLength / 8
        if (n == 1) {
            return new DataView(ab).getFloat64(0, true);
        }
        var arr = new Array(n);
        var dview = new DataView(ab);
        for (var i = 0; i < n; i++) {
            arr[i] = dview.getFloat64(i * 8, true);
        }
        return arr
    }

    var res = {}
    res.header = {
        abif: RTC(rawdata.slice(0, 4)),
        version: SInt16(rawdata.slice(4, 6)),
        name: RTC(rawdata.slice(6, 10)),
        number: SInt32(rawdata.slice(10, 14)),
        elementtype: SInt16(rawdata.slice(14, 16)),
        elementsize: SInt16(rawdata.slice(16, 18)),
        numelements: SInt32(rawdata.slice(18, 22)),
        dataoffset: SInt32(rawdata.slice(26, 30)),
        datahandle: SInt32(rawdata.slice(30, 34))
    }
    res.directory = []
    res.data = {}

    for (var i = 0; i < res.header.numelements; i++) {
        var deb = i * res.header.elementsize + res.header.dataoffset
        var direntry = rawdata.slice(deb, deb + res.header.elementsize)
        var name = RTC(direntry.slice(0, 4))
        var tagnumber = SInt32(direntry.slice(4, 8));
        var elementtype = name == "PCON" ? 1 : SInt16(direntry.slice(8, 10))
        var elementsize = SInt16(direntry.slice(10, 12));
        var datasize = SInt32(direntry.slice(16, 20));
        var numelements = SInt32(direntry.slice(12, 16));
        var dataoffset = 0;
        var data;
        if (datasize <= 4) {
            dataoffset = deb + 20
        } else {
            dataoffset = SInt32(direntry.slice(20, 24));
            data = rawdata.slice(dataoffset, dataoffset + numelements * elementsize)
        }
        data = rawdata.slice(dataoffset, dataoffset + numelements * elementsize)

        res.directory.push({
            name: name,
            tagnumber: tagnumber,
            elementtype: elementtype,
            elementsize: elementsize,
            numelements: numelements,
            datasize: datasize,
            dataoffset: dataoffset
        })

        if (elementtype == 1) {
            data = UInt8(data, numelements)
        } else if (elementtype == 2) {
            data = RTC(data)
        } else if (elementtype == 3) {
            data = UInt16(data, numelements)
        } else if (elementtype == 4) {
            data = SInt16(data, numelements)
        } else if (elementtype == 5) {
            data = SInt32(data, numelements)
        } else if (elementtype == 7) {
            data = f32(data, numelements)
        } else if (elementtype == 8) {
            data = f64(data, numelements)
        } else if (elementtype == 10) {
            data = {
                year: SInt16(data.slice(0, 2)),
                month: UInt8(data.slice(2, 3)),
                day: UInt8(data.slice(3, 4))
            }
        } else if (elementtype == 11) {
            data = {
                hour: UInt8(data.slice(0, 1)),
                minute: UInt8(data.slice(1, 2)),
                second: UInt8(data.slice(2, 3)),
                hsecond: UInt8(data.slice(3, 4)),
            }
        } else if (elementtype == 18) {
            data = RTC(data.slice(1))
        } else if (elementtype == 19) {
            data = RTC(data.slice(0, -1))
        } else if (elementtype >= 1024) {
            data = new Uint8Array(data)
        } else {
            throw "unsupported type found in file"
        }
        res.data[name + "_" + tagnumber] = data
    }

    return res
}
