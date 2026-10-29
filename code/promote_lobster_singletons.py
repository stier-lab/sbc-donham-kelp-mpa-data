"""Promote only lobster biomass after repairing singleton resampling."""
import argparse
import json
from pathlib import Path
from promote_monthly_kelp import FILES, md5, prepare


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--target', type=Path, action='append', required=True)
    p.add_argument('--record', type=Path, required=True)
    p.add_argument('--apply', action='store_true')
    p.add_argument('--urchins', action='store_true')
    a = p.parse_args()
    products, results, details = {}, {}, {}
    for name in FILES:
        old = (a.baseline/name).read_bytes()
        products[name], results[name], details[name] = prepare(old,(a.candidate/name).read_bytes(),
            'raw' in name,taxon_name=('Mesocentrotus franciscanus','Strongylocentrotus purpuratus') if a.urchins else 'Panulirus interruptus',
            sources=('PISCO','LTER') if a.urchins else ('PISCO',),responses=('Bio',))
        for target in a.target:
            if target.resolve() == a.baseline.resolve() or (target/name).read_bytes() != old:
                raise ValueError('Destination does not match separate baseline: '+str(target/name))
    report = {'release':'lobster_singleton_sampling_2026_10_18','applied':a.apply,
        'scope':{'programme':'PISCO','taxon':'Panulirus interruptus','response':'Bio'},
        'sampling':'Sample vector positions, not singleton numeric values; only singleton cached biomass replaced',
        'a':0.001352821,'b':2.913963113,'length_unit':'mm','mass_unit':'g',
        'files':results}
    if a.urchins:
        report['release'] = 'urchin_singleton_sampling_2026_10_18'
        report['scope'] = {'programmes':['PISCO','LTER'],
            'taxa':['Mesocentrotus franciscanus','Strongylocentrotus purpuratus'],'response':'Bio'}
        report['sampling'] = 'Sample vector positions; replace only singleton cached biomass totals, preserving non-singleton draws'
        for field in ('a','b','length_unit','mass_unit'):
            del report[field]
    print(json.dumps(report,indent=2))
    if a.apply:
        for target in a.target:
            for name, data in products.items():
                (target/name).write_bytes(data)
            marker = 'urchin_bootstrap_correction.json' if a.urchins else 'lobster_bootstrap_correction.json'
            (target/marker).write_text(json.dumps(report,indent=2)+'\n')
            (target/'checksums.md5').write_text(''.join(md5(f.read_bytes())+' '+f.name+'\n' for f in sorted(target.glob('*.csv'))))
        a.record.write_text(json.dumps({**report,'changed_records':details},indent=2)+'\n')


if __name__ == '__main__':
    main()
