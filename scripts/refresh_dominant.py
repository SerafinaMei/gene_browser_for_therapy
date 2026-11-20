"""Recompute a conservative dominant P/LP allele subset from the weekly ClinVar VCV XML.

Counts require a current, single-gene SimpleAllele, a P/LP variant-condition
classification with assertion criteria, and resolved dominant inheritance.
Disease mappings are condition-specific, never inherited from a gene-wide label.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,re,subprocess,tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'dist/data'
CACHE=ROOT/'data_sources/dominant-current.json.gz'
XML_URL='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/ClinVarVCVRelease_00-latest_weekly.xml.gz'
CONDITION_URL='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id'
PLP={'Pathogenic','Likely pathogenic','Pathogenic/Likely pathogenic'}
REVIEWS={'criteria provided, single submitter','criteria provided, multiple submitters, no conflicts','reviewed by expert panel','practice guideline'}
DOM={'AD','XLD','SD'}
HEADERS=['AlleleID','VariationID','Variant','Status','RCV','ConditionIDs','Condition','Classification','ReviewStatus','InheritanceEvidence']

def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,separators=(',',':'),allow_nan=False))

def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def fetch(url,path):
    subprocess.run(['curl','--fail','--silent','--show-error','--retry','3','--connect-timeout','20','--max-time','5400',url,'-o',str(path)],check=True)

def modes(text):
    text=(text or '').lower()
    found=set()
    for label,code in [('autosomal dominant','AD'),('autosomal recessive','AR'),('x-linked dominant','XLD'),('x-linked recessive','XLR'),('semidominant','SD'),('mitochondrial','MT'),('y-linked','YL')]:
        if label in text:found.add(code)
    if text.strip().upper() in {'AD','AR','XLD','XLR','SD','MT','YL','XL'}:found.add(text.strip().upper())
    if 'x-linked' in text and not found:found.add('XL')
    return found

def key(db,value):
    db=(db or '').upper();value=value or ''
    if db in {'OMIM','MIM'}:return 'OMIM:'+value.removeprefix('OMIM:')
    if db=='MONDO':return 'MONDO:'+value.removeprefix('MONDO:')
    if db=='MEDGEN':return 'MedGen:'+value
    return None

def identifiers(element):
    if element is None:return set()
    return {k for x in element.iter('XRef') if (k:=key(x.get('DB'),x.get('ID')))}

def load_context(condition_file):
    genes={}
    for p in DATA.glob('genes-*.json'):genes.update(json.loads(p.read_text()))
    by_id={str(g['entrez_id']):s for s,g in genes.items() if g.get('entrez_id')}
    lookup=defaultdict(set)
    for symbol,g in genes.items():
        for r in g['diseases']:
            if r.get('CLASSIFICATION') in {'Definitive','Strong','Moderate'}:
                lookup[(symbol,r['DISEASE ID (MONDO)'])].update(modes(r.get('MOI')))
    for p in DATA.glob('omim-*.json'):
        for symbol,wrapper in json.loads(p.read_text()).items():
            for r in wrapper['omim_relationships']['records']:
                if r['phenotype_mim'] and r['mapping_key']==3:
                    lookup[(symbol,'OMIM:'+r['phenotype_mim'])].update(modes(r['inheritance']))
    aliases=defaultdict(set)
    with condition_file.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not {'#GeneID','ConceptID','DiseaseMIM','SourceName','SourceID'}.issubset(reader.fieldnames or []):raise ValueError('ClinVar condition mapping schema changed')
        for r in reader:
            symbol=by_id.get(r['#GeneID'])
            if not symbol:continue
            ids={x for x in [key('MedGen',r['ConceptID']) if r['ConceptID'] else None,key(r['SourceName'],r['SourceID']) if r['SourceID'] else None,key('OMIM',r['DiseaseMIM']) if re.fullmatch(r'\d{6}',r['DiseaseMIM']) else None] if x}
            for identifier in ids:aliases[(symbol,identifier)].update(ids)
    return genes,by_id,lookup,aliases

def expand(symbol,ids,aliases):
    return set(ids).union(*(aliases.get((symbol,i),set()) for i in ids))

def evaluate_record(record,by_id,lookup,aliases):
    if record.get('RecordType')!='classified' or record.findtext('RecordStatus')!='current':return None
    classified=record.find('ClassifiedRecord');allele=classified.find('SimpleAllele') if classified is not None else None
    if allele is None:return None # Haplotype/genotype and included-only records are not individual classified alleles.
    aggregate=classified.find('Classifications/GermlineClassification')
    if aggregate is None or aggregate.findtext('Description') not in PLP:return None
    if aggregate.findtext('ReviewStatus') not in REVIEWS:return None
    gene_nodes=allele.findall('GeneList/Gene')
    if len(gene_nodes)!=1:return None # Avoid assigning overlapping multigene structural alleles to every gene.
    symbol=by_id.get(gene_nodes[0].get('GeneID'))
    if not symbol:return None
    allele_id=allele.get('AlleleID')
    if not allele_id:raise ValueError('Missing AlleleID')
    # Join submitter inheritance to the same condition via explicit identifiers or ClinVar's TraitMapping.
    mapped=defaultdict(set)
    for mapping in classified.findall('TraitMappingList/TraitMapping'):
        for medgen in mapping.findall('MedGen'):
            if medgen.get('CUI'):mapped[mapping.get('ClinicalAssertionID')].add('MedGen:'+medgen.get('CUI'))
    direct=[]
    for assertion in classified.findall('ClinicalAssertionList/ClinicalAssertion'):
        if assertion.get('ContributesToAggregateClassification')=='false':continue
        if assertion.findtext('RecordStatus') not in ('current',None):continue
        if assertion.findtext('Classification/GermlineClassification') not in PLP:continue
        if assertion.findtext('Classification/ReviewStatus') not in REVIEWS:continue
        inheritance=set()
        for a in assertion.findall('AttributeSet/Attribute'):
            if a.get('Type')=='ModeOfInheritance':inheritance.update(modes(a.text))
        if inheritance:
            ids=identifiers(assertion.find('TraitSet'))|mapped[assertion.get('ID')]
            direct.append((expand(symbol,ids,aliases),inheritance))
    rows=[]
    for rcv in classified.findall('RCVList/RCVAccession'):
        classification=rcv.find('RCVClassifications/GermlineClassification')
        if classification is None or classification.findtext('Description') not in PLP:continue
        if classification.findtext('ReviewStatus') not in REVIEWS:continue
        conditions=rcv.findall('ClassifiedConditionList/ClassifiedCondition')
        condition_ids={k for c in conditions if (k:=key(c.get('DB'),c.get('ID')))}
        # Resolve every condition separately; a known AD condition cannot resolve
        # another unnamed or unmapped condition in the same RCV record.
        evidence=set(); explicit_used=False; missing=not bool(conditions)
        for condition in conditions:
            identifier=key(condition.get('DB'),condition.get('ID'))
            expanded=expand(symbol,{identifier} if identifier else set(),aliases)
            explicit=set().union(*(m for ids,m in direct if ids & expanded))
            inferred=set().union(*(lookup.get((symbol,i),set()) for i in expanded))
            resolved=explicit or inferred
            missing=missing or not bool(resolved)
            explicit_used=explicit_used or bool(explicit)
            evidence.update(resolved)
        # Mixed modes and unspecified X-linked inheritance are unresolved.
        if not missing and evidence and evidence<=DOM:
            status='dominant'
            basis=('ClinVar variant-condition inheritance and/or condition mapping: ' if explicit_used else 'Condition-specific inheritance mapping: ')+','.join(sorted(evidence))
        elif not missing and evidence and evidence<={'AR','XLR','MT','YL'}:
            continue
        else:
            status='inheritance_unresolved';basis='Mixed or unavailable inheritance'
        rows.append([allele_id,record.get('VariationID'),record.get('VariationName'),status,rcv.get('Accession'),';'.join(sorted(condition_ids)),'; '.join(c.text or '' for c in conditions),classification.findtext('Description'),classification.findtext('ReviewStatus'),basis])
    return (symbol,rows) if rows else None

def parse_release(xml_path,context,minimum_records=1000000):
    genes,by_id,lookup,aliases=context
    result=defaultdict(list);records=0
    opener=gzip.open if xml_path.suffix=='.gz' else open
    with opener(xml_path,'rb') as f:
        events=ET.iterparse(f,events=('start','end'))
        _,root=next(events)
        date=root.get('ReleaseDate')
        if root.tag!='ClinVarVariationRelease' or not date:raise ValueError('Unexpected ClinVar release schema')
        for event,node in events:
            if event=='end' and node.tag=='VariationArchive':
                item=evaluate_record(node,by_id,lookup,aliases)
                if item:result[item[0]].extend(item[1])
                records+=1;root.clear()
                if records%250000==0:
                    position=f.fileobj.tell() if isinstance(f,gzip.GzipFile) else f.tell()
                    print(f'Parsed {records:,} ClinVar records ({position/xml_path.stat().st_size:.0%} of input)',flush=True)
    if records<minimum_records:raise ValueError(f'Incomplete ClinVar release: only {records} records')
    if not result:raise ValueError('No eligible P/LP records: check ClinVar schema')
    source=dict(date=date,url=XML_URL,sha256=sha(xml_path),records=records,checked_at=datetime.now(timezone.utc).isoformat())
    source['clinvar_condition_map_sha256']=sha(ROOT/'data_sources/clinvar-gene-conditions.tsv') if (ROOT/'data_sources/clinvar-gene-conditions.tsv').exists() else None
    source['omim_date']=json.loads((ROOT/'data_sources/omim-gene-phenotypes.json').read_text())['source']['date']
    source['clingen_date']=next(s['date'] for s in json.loads((DATA/'manifest.json').read_text())['sources'] if s['id']=='validity')
    return dict(source=source,headers=HEADERS,genes=dict(result))

def build(data):
    source=data['source'];counts={}
    for path in DATA.glob('genes-*.json'):
        letter=path.stem.removeprefix('genes-');shard={};columns={}
        for symbol,gene in json.loads(path.read_text()).items():
            records=data['genes'].get(symbol,[])
            dominant={r[0] for r in records if r[3]=='dominant'}
            unresolved={r[0] for r in records if r[3]=='inheritance_unresolved'}-dominant
            assessed=bool(gene.get('entrez_id'))
            count=dict(count=len(dominant) if assessed else None,unresolved=len(unresolved) if assessed else None,date=source['date'],assessed=assessed)
            if records:
                columns[symbol]=[list(col) for col in zip(*records)]
                count['export']=dict(path=f'downloads/dominant/{letter}.json.gz',format='columnar-gzip',record_key=symbol)
            shard[symbol]={'dominant_current':count};counts[symbol]=count
        save(DATA/f'dominant-{letter}.json',shard)
        path=ROOT/'dist/downloads/dominant'/f'{letter}.json.gz';path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(gzip.compress(json.dumps(dict(headers=HEADERS,genes=columns),separators=(',',':')).encode(),mtime=0))
    path=DATA/'index.json';index=json.loads(path.read_text());ci=index['fields'].index('dominant_plp');si=index['fields'].index('symbol')
    for row in index['rows']:row[ci]=counts[row[si]]['count']
    save(path,index)
    manifest=json.loads((DATA/'manifest.json').read_text())
    note='Distinct single-gene ClinVar alleles with P/LP germline variant and variant-condition classifications and assertion criteria at both levels. Dominant inheritance comes from matching ClinVar submissions or exact condition-ID mappings to ClinGen/OMIM. Mixed or missing inheritance is excluded and reported separately. Overlapping multigene alleles and haplotype/genotype records are excluded. OMIM inheritance snapshot: '+source['omim_date']+'. ClinGen inheritance snapshot: '+source['clingen_date']+'. Genes without an NCBI Gene ID are unavailable, not zero.'
    manifest['sources']=[s for s in manifest['sources'] if s['id']!='dominant_current']+[dict(id='dominant_current',title='Current dominant P/LP variants',url=source['url'],date=source['date'],path=str(CACHE.relative_to(ROOT)),sha256=sha(CACHE),note=note)]
    save(DATA/'manifest.json',manifest)
    print(f"Dominant subset: {sum(c['count'] or 0 for c in counts.values()):,} alleles across {sum(bool(c['count']) for c in counts.values()):,} genes",flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--fetch',action='store_true');parser.add_argument('--xml',type=Path);parser.add_argument('--work-dir',type=Path,default=Path(tempfile.gettempdir())/'dnd-weekly');args=parser.parse_args()
    mapping=ROOT/'data_sources/clinvar-gene-conditions.tsv'
    if args.fetch:
        args.work_dir.mkdir(parents=True,exist_ok=True)
        xml=args.work_dir/'clinvar-weekly.xml.gz';checksum=args.work_dir/'clinvar.md5'
        fetch(XML_URL+'.md5',checksum);expected=checksum.read_text().split()[0]
        fetch(XML_URL,xml)
        with xml.open('rb') as f:actual=hashlib.file_digest(f,'md5').hexdigest()
        if actual!=expected:raise ValueError('ClinVar checksum changed or download incomplete')
        fetch(CONDITION_URL,mapping);args.xml=xml
    if args.xml:
        data=parse_release(args.xml,load_context(mapping))
        CACHE.write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0))
    else:
        with gzip.open(CACHE,'rt') as f:data=json.load(f)
    build(data)
