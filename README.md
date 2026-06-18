# abif.js

A javascript function to parse Applied Biosystems ABIF files

## usage

```
const {
    header,     // file metadata
    directory,  // fields metadata
    data        // data fields
} = parse_abif(await ab1file.arrayBuffer())

let fasta = `>${data.SMPL_1}
${data.PBAS_1.match(/.{1,70}/g).join("\n")}`
let fastq = `@${data.SMPL_1}
${data.PBAS_1}
+
${data.PCON_1.map(x => String.fromCharCode(x + 33)).join("")}`
```
