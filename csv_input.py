"""Native diffraction-table input validation; no private framework dependencies."""
import csv,io,math,pathlib
import numpy as np

ALIASES={'h':'h','k':'k','l':'l','2theta':'angle','two_theta_deg':'angle','d-value':'d','d_value':'d','mult':'mult','intensity':'intensity'}

def load_diffraction_table(path,require_intensity=True,return_delimiter=False):
    text=pathlib.Path(path).read_text(encoding='utf-8-sig')
    lines=[line for line in text.splitlines() if line.strip()]
    if not lines:raise ValueError('Empty diffraction table.')
    delimiter=';' if ';' in lines[0] else ','
    reader=csv.DictReader(io.StringIO('\n'.join(lines)),delimiter=delimiter)
    names={};seen=set()
    for name in reader.fieldnames or []:
        canonical=ALIASES.get(name.strip().lower())
        if canonical:
            if canonical in seen:raise ValueError('Duplicate column: '+canonical)
            seen.add(canonical);names[canonical]=name
    needed={'h','k','l','angle'}|({'intensity'} if require_intensity else set())
    missing=needed-set(names)
    if missing:raise ValueError('Missing columns: '+', '.join(sorted(missing)))
    rows=[]
    for number,row in enumerate(reader,2):
        if None in row:raise ValueError(f'Row {number}: extra fields.')
        def number_value(key,default=0.):
            if key not in names:return default
            value=row.get(names[key])
            if value is None or not value.strip():raise ValueError(f'Row {number}: missing {key}.')
            if delimiter==';':value=value.replace(',','.')
            try:value=float(value)
            except ValueError as error:raise ValueError(f'Row {number}: nonnumeric {key}.') from error
            if not math.isfinite(value):raise ValueError(f'Row {number}: nonfinite {key}.')
            return value
        h,k,l,angle=[number_value(key) for key in ('h','k','l','angle')]
        if not all(v.is_integer() for v in (h,k,l)):raise ValueError(f'Row {number}: Miller indices must be integers.')
        if h==k==l==0:raise ValueError(f'Row {number}: 000 is not a diffraction reflection.')
        if not 0<angle<180:raise ValueError(f'Row {number}: require 0 < 2theta < 180 degrees.')
        rows.append([h,k,l,angle,number_value('d'),number_value('mult'),number_value('intensity')])
    if not rows:raise ValueError('No diffraction rows.')
    data=np.asarray(rows,dtype=float).reshape(-1,7)
    return (data,delimiter) if return_delimiter else data
