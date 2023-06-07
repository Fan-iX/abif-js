# abif.js

A javascript script for parsing Applied Biosystems ABIF file

## usage

```
parse_abif(arraybuffer: rawdata)
```

Returns a object:

```
{
    header: {...}, // file metadata
    directory: [...], // fields metadata
    data: {...} // data fields
}
```
