"""Decode a bounded NASA/PDS ODF sample using its own PDS4 label.

Only inventory, layout, time tags, stations and type flags are interpreted.
This is not an orbit solution, calibrated residual series, or Q2 estimate.
"""
from pathlib import Path
import xml.etree.ElementTree as E
import json, hashlib, csv
from collections import Counter
import numpy as np
R=Path(__file__).resolve().parent
NS={'p':'http://pds.nasa.gov/pds4/pds/v1'}
def txt(e,path):return e.findtext(path,namespaces=NS)
def main():
    p=R/'raw/cassini_sagr_sample.dat';data=p.read_bytes()
    root=E.parse(R/'raw/cassini_sagr_sample.xml').getroot()
    assert len(data)==int(txt(root,'.//p:File/p:file_size'))
    assert hashlib.md5(data).hexdigest()==txt(root,'.//p:File/p:md5_checksum')
    tables=[];fields={};orbit=None;headers=[]
    for tb in root.findall('.//p:Table_Binary',NS):
        name=txt(tb,'p:name');offset=int(txt(tb,'p:offset'));count=int(txt(tb,'p:records'));length=int(txt(tb,'p:Record_Binary/p:record_length'))
        assert length==36 and offset+count*length<=len(data)
        tables.append(dict(name=name,offset=offset,records=count,record_length=length))
        if 'Header' in name or 'End-of-File' in name:
            key=int.from_bytes(data[offset:offset+4],'big',signed=True)
            headers.append(dict(name=name,primary_key=key))
        if name=='ODF File Label Group Data':
            for f in tb.findall('p:Record_Binary/p:Field_Binary',NS):
                pos=offset+int(txt(f,'p:field_location'))-1;n=int(txt(f,'p:field_length'));dt=txt(f,'p:data_type')
                fields[txt(f,'p:name')]=data[pos:pos+n].decode().strip() if dt=='ASCII_String' else int.from_bytes(data[pos:pos+n],'big')
        if name=='ODF Orbit Data Group Data':orbit=tb
    assert headers[-1]['primary_key']==-1
    offset=int(txt(orbit,'p:offset'));count=int(txt(orbit,'p:records'))
    words=np.frombuffer(data[offset:offset+count*36],dtype='>u4').reshape(count,9).astype(np.uint64)
    bit_defs={}
    for field in orbit.findall('p:Record_Binary/p:Field_Binary',NS):
        loc=int(txt(field,'p:field_location'))-1;n=int(txt(field,'p:field_length'))
        if n!=4:continue
        for bit in field.findall('p:Packed_Data_Fields/p:Field_Bit',NS):
            start=int(txt(bit,'p:start_bit_location'));stop=int(txt(bit,'p:stop_bit_location'))
            bit_defs[txt(bit,'p:name')]=(words[:,loc//4]>>(32-stop))&((1<<(stop-start+1))-1)
    ms=bit_defs['Record Time Tag, fractional part'];assert np.all(ms<1000)
    stamps=words[:,0].astype(float)+ms/1000
    names=['Receiving Station ID','Transmitting Station ID','Data Type ID','Data Validity Indicator']
    inventories={}
    for name in ['sagr','sce']:
        rows=list(csv.DictReader((R/f'raw/cassini_{name}_inventory.csv').open()))
        inventories[name]=dict(rows=len(rows),columns=list(rows[0]))
    result=dict(scope='Bounded archive accessibility and binary decoding, not Q2 fitting',
        source_filename=txt(root,'.//p:File/p:file_name'),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),
        label_md5_verified=True,file_label=fields,tables=tables,header_keys=headers,
        orbit_records=count,time_tags_nondecreasing=bool(np.all(np.diff(stamps)>=0)),
        time_seconds_from_file_reference=[float(stamps.min()),float(stamps.max())],
        categorical_counts={n:dict(Counter(map(int,bit_defs[n]))) for n in names},inventories=inventories,
        remaining_bytes_after_logical_eof=len(data)-tables[-1]['offset']-36)
    np.savetxt(R/'results/cassini_odf_decoded_index.csv',np.column_stack([stamps]+[bit_defs[n] for n in names]),delimiter=',',
        header='seconds_from_file_reference,receiving_station,transmitting_station,data_type,validity_flag',comments='',fmt=['%.3f','%d','%d','%d','%d'])
    (R/'results/cassini_odf_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
