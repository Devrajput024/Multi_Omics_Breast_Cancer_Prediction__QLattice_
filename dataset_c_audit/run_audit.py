from pathlib import Path
import csv, json, math, hashlib, collections, statistics, re, zipfile, html

ROOT = Path.cwd()
EXTRACT = ROOT / 'dataset_c_extracted'
OUT = ROOT / 'dataset_c_audit'
ZIP = Path(r'C:\Users\asus\Downloads\archive (8).zip')
FILES = ['data.csv', 'brca_data_w_subtypes.csv']

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        r=csv.reader(f); header=next(r); rows=list(r)
    return header, rows

def percentile(xs, p):
    if not xs: return ''
    ys=sorted(xs); pos=(len(ys)-1)*p; lo=int(pos); hi=min(lo+1,len(ys)-1); frac=pos-lo
    return ys[lo]*(1-frac)+ys[hi]*frac

def numeric_column(vals):
    out=[]
    for v in vals:
        if v.strip()=='': continue
        try: out.append(float(v))
        except ValueError: return None
    return out

def safe_num(x):
    if x=='': return ''
    if not math.isfinite(x): return ''
    return round(x,10)

def pearson(x,y):
    n=len(x)
    if n<2: return 0.0
    mx=sum(x)/n; my=sum(y)/n
    xx=sum((a-mx)**2 for a in x); yy=sum((b-my)**2 for b in y)
    if xx==0 or yy==0: return 0.0
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(xx*yy)

def esc(s): return html.escape(str(s),quote=True)

def svg_write(name, body, width=1000, height=600):
    svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><style>text{{font-family:Arial,sans-serif;fill:#222}} .small{{font-size:12px}} .label{{font-size:14px}} .title{{font-size:20px;font-weight:bold}}</style>{body}</svg>'
    (OUT/'figures'/name).write_text(svg,encoding='utf-8')

def textline(s): return s+'\n'

# Load only the two tabular files contained in the selected Dataset C archive.
data={fn:read_csv(EXTRACT/fn) for fn in FILES}
profiles={}
with (OUT/'dataset_structure.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f)
    w.writerow(['dataset_file','column_name','inferred_type','role','unique_nonmissing','missing_count','missing_percent','constant','near_constant_95pct','mean','std_sample','min','q25','median','q75','max','first5_values_sha256','last5_values_sha256'])
    for fn,(header,rows) in data.items():
        n=len(rows); p=len(header); colstats=[]; colvecs={}
        for j,col in enumerate(header):
            vals=[r[j] if j<len(r) else '' for r in rows]
            non=[v for v in vals if v.strip()!='']; miss=n-len(non)
            nums=numeric_column(vals)
            kind='numeric' if nums is not None else 'categorical'
            uniques=len(set(non)); counts=collections.Counter(non)
            modal=max(counts.values(),default=0); constant=(uniques<=1); near=(len(non)>0 and modal/len(non)>=.95)
            if col=='vital.status': role='target used in notebook; clinical meaning undocumented'
            elif col in ['PR.Status','ER.Status','HER2.Final.Status','histological.type']: role='clinical/phenotype candidate; exact definitions undocumented'
            else: role='feature; biological modality not established by data dictionary'
            first_hash=hashlib.sha256('\x1f'.join(vals[:5]).encode()).hexdigest()
            last_hash=hashlib.sha256('\x1f'.join(vals[-5:]).encode()).hexdigest()
            if nums is not None:
                mean=sum(nums)/len(nums) if nums else 0
                sd=statistics.stdev(nums) if len(nums)>1 else 0
                q1=percentile(nums,.25); med=percentile(nums,.5); q3=percentile(nums,.75)
                lo=percentile(nums,.25); hi=percentile(nums,.75); iqr=hi-lo
                outliers=sum(v < lo-1.5*iqr or v > hi+1.5*iqr for v in nums)
                colvecs[col]=nums
            else:
                mean=sd=q1=med=q3=0
            w.writerow([fn,col,kind,role,uniques,miss,round(100*miss/n,8),str(constant).lower(),str(near).lower(),safe_num(mean),safe_num(sd),safe_num(min(nums) if nums else 0) if nums is not None else '',safe_num(q1),safe_num(med),safe_num(q3),safe_num(max(nums) if nums else 0) if nums is not None else '',first_hash,last_hash])
            colstats.append({'name':col,'missing':miss,'unique':uniques,'constant':constant,'near':near,'type':kind,'mode_share':modal/len(non) if non else 0,'outliers':outliers if nums is not None else 0})
        profiles[fn]={'n':n,'p':p,'header':header,'rows':rows,'stats':colstats,'vectors':colvecs}

# archive inventory, file types, and sizes
with zipfile.ZipFile(ZIP) as z:
    infos=z.infolist()
    inv=[]
    for x in infos:
        suffix=Path(x.filename).suffix.lower() or '[no extension]'
        inv.append((x.filename,x.file_size,suffix,'directory' if x.is_dir() else 'file'))
with (OUT/'dataset_inventory.txt').open('w',encoding='utf-8') as f:
    f.write(f'Dataset C archive: {ZIP}\nArchive size: {ZIP.stat().st_size:,} bytes\nExtracted to: {EXTRACT}\n\nArchive entries ({len(inv)}):\n')
    for path,size,fmt,typ in inv: f.write(f'- {path} | {typ} | {fmt} | {size:,} bytes\n')
    f.write('\nArchive file summary:\n')
    f.write('- CSV: 2 (brca_data_w_subtypes.csv, data.csv)\n- Jupyter notebook: 1\n- Markdown documentation: 1\n- QLattice model artifacts: 2\n- XLSX/TXT/PDF/JSON files: 0\n- Explicit directory entries: 0 (all files at archive root)\n')

# per-column missingness
with (OUT/'missing_values.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(['dataset_file','column_name','missing_count','missing_percent_of_column','missing_percent_of_all_cells'])
    for fn,pr in profiles.items():
        for s in pr['stats']:
            w.writerow([fn,s['name'],s['missing'],round(100*s['missing']/pr['n'],8),round(100*s['missing']/(pr['n']*pr['p']),10)])

# Target analysis
br=profiles['brca_data_w_subtypes.csv']; dc=profiles['data.csv']
targets=['vital.status','PR.Status','ER.Status','HER2.Final.Status','histological.type']
with (OUT/'target_analysis.txt').open('w',encoding='utf-8') as f:
    f.write('TARGET VARIABLE ANALYSIS (computed from Dataset C files only)\n\n')
    f.write('Documentation finding: README describes potential subtype/prognosis objectives generically but does not define actual labels. Notebook code explicitly sets output_name="vital.status" with kind="classification". Both saved .model artifacts also encode output name vital.status. No codebook defines what 0 or 1 means, time point, censoring, or event mapping. Do not translate codes to alive/dead or survival without source documentation.\n\n')
    for fn,pr in profiles.items():
        f.write(f'{fn} ({pr["n"]} rows):\n')
        for target in targets:
            if target not in pr['header']: continue
            j=pr['header'].index(target); vals=[r[j] if r[j].strip() else '<MISSING>' for r in pr['rows']]; cnt=collections.Counter(vals)
            f.write(f'  {target}: missing={cnt.get("<MISSING>",0)}; observed classes/categories:\n')
            for k,v in sorted(cnt.items()): f.write(f'    {k}: {v} ({100*v/pr["n"]:.4f}%)\n')
            if target=='vital.status':
                observed=[(k,v) for k,v in cnt.items() if k!='<MISSING>']; maj=max((v for k,v in observed),default=0); mi=min((v for k,v in observed),default=0)
                f.write(f'    majority/minority ratio={maj/mi:.4f}:1; classification imbalance is present. Classification is notebook-intended, but target semantics remain unverified. Regression is unsupported by this binary-coded variable as supplied.\n')
        f.write('\n')
    f.write('Potential leakage: vital.status is present as an input column in both files; it must be removed from predictors when it is the target. The subtype-enriched file adds receptor-status and histology labels; if one of those is selected as outcome, related clinical markers and post-outcome measurements may leak information depending on timing. The archive provides no measurement dates or prediction-time definition. Feature/target chronology therefore cannot be cleared.\n')

# Archive documentation, notebook, saved model facts
readme=next(EXTRACT.glob('*readme*')).read_text(encoding='utf-8',errors='replace')
nb=json.loads(next(EXTRACT.glob('*.ipynb')).read_text(encoding='utf-8'))
code='\n'.join(''.join(c.get('source',[])) for c in nb.get('cells',[]) if c.get('cell_type')=='code')
model_facts=[]
for mp in EXTRACT.glob('*.model'):
    try:
        obj=json.loads(mp.read_text(encoding='utf-8'))
        if isinstance(obj,dict):
            names=[]; outnames=[]
            if 'nodes' in obj:
                names=[x.get('name','') for x in obj['nodes'] if x.get('name')]
                outnames=[x.get('name','') for x in obj['nodes'] if str(x.get('spec','')).startswith('out:')]
            elif 'program' in obj:
                names=[x for x in obj.get('names',[]) if x]
                outnames=[x for x,fnm in zip(obj.get('names',[]),obj.get('fnames',[])) if str(fnm).startswith('out:')]
            model_facts.append((mp.name, outnames, names))
    except Exception: pass
# Feature/mode mapping report
with (OUT/'feature_category_analysis.txt').open('w',encoding='utf-8') as f:
    f.write('FEATURE CATEGORY AND MODALITY INVESTIGATION\n\n')
    f.write('All feature family labels below are name-based observations, not confirmed biological modality assignments. The archive has no feature dictionary, assay specification, or external sample manifest. The README gives broad examples of expression, mutation/CNV, and clinical data, but does not map these prefixes or provide layer definitions.\n\n')
    prefixes=collections.Counter(c.split('_',1)[0] for c in dc['header'] if '_' in c)
    examples=collections.defaultdict(list)
    for c in dc['header']:
        if '_' in c:
            prefix=c.split('_',1)[0]
            if len(examples[prefix])<8: examples[prefix].append(c)
    for key in ['rs','mu','cn','pp']:
        f.write(f'{key}_ family: {prefixes.get(key,0)} columns; examples: {", ".join(examples.get(key,[]))}. Numeric values in both tables. Biological meaning not verified from documentation.\n')
    f.write('\nCautious interpretation: `mu_` names and saved model feature names are consistent with mutation-labelled variables; `cn_` names are consistent with copy-number-labelled variables; `pp_` names look protein/phosphoprotein-like; `rs_` cannot be safely called SNP or transcriptomic data from its prefix alone. No standalone demographic columns were found. `PR.Status`, `ER.Status`, `HER2.Final.Status`, `histological.type`, and `vital.status` are clinical/phenotype-like labels by name, but their codebooks are absent.\n\n')
    f.write('Files store all four numeric feature families in one wide matrix. No separate omics files exist in the archive. All 1,936 non-target columns in data.csv and the corresponding shared columns in the enriched table are numeric; the four extra clinical fields are categorical.\n')

# sample matching and feature alignment
h1,r1=br['header'],br['rows']; h2,r2=dc['header'],dc['rows']
shared=[c for c in h2 if c in set(h1)]
idx1=[h1.index(c) for c in shared]; idx2=[h2.index(c) for c in shared]
exact_rows=sum(all(a[i]==b[j] for i,j in zip(idx1,idx2)) for a,b in zip(r1,r2)) if len(r1)==len(r2) else 0
id_cols=[c for c in set(h1)|set(h2) if re.search(r'(^|[._])(patient|sample|case|submitter|barcode|id)([._]|$)',c,re.I)]
with (OUT/'sample_overlap_report.txt').open('w',encoding='utf-8') as f:
    f.write('SAMPLE OVERLAP AND INTEGRATION READINESS\n\n')
    f.write(f'Tabular files: {len(FILES)}. Rows: data.csv={len(r2)}, brca_data_w_subtypes.csv={len(r1)}. Columns shared: {len(shared)}. Additional columns in brca_data_w_subtypes.csv: {", ".join(c for c in h1 if c not in set(h2))}.\n')
    f.write(f'Exact row-wise matches across shared columns at same position: {exact_rows}/{min(len(r1),len(r2))}. All shared headers are in the same order: {h1[:len(h2)]==h2}.\n')
    f.write(f'Explicit patient/sample ID columns detected by header name: {len(id_cols)}.\n')
    f.write('No patient/sample IDs exist in the tables, so matched/unmatched patient counts and duplicate IDs are not computable. Exact row order alignment is not equivalent to verified patient ID matching. The added labels align by row position for these two supplied tables; do not use row order to join with external cohorts.\n')
    f.write('Feature blocks can be used together in a single aligned matrix as stored, but whether they are measurements from the same patients cannot be proven independently. No automatic merge was performed.\n')

# Quality report and data-level EDA
with (OUT/'data_quality_report.txt').open('w',encoding='utf-8') as f:
    f.write('DATA QUALITY AND STATISTICAL EDA (Dataset C only)\n\n')
    for fn,pr in profiles.items():
        cells=pr['n']*pr['p']; miss=sum(s['missing'] for s in pr['stats']); dup_rows=pr['n']-len({tuple(r) for r in pr['rows']})
        colgroups=collections.defaultdict(list)
        for j,c in enumerate(pr['header']): colgroups[tuple(r[j] for r in pr['rows'])].append(c)
        dups=[v for v in colgroups.values() if len(v)>1]
        f.write(f'{fn}: {pr["n"]} rows x {pr["p"]} columns; total missing cells={miss}/{cells} ({100*miss/cells:.6f}%); duplicate rows={dup_rows}; identical-value column groups={len(dups)} (extra redundant columns={sum(len(x)-1 for x in dups)}); constant columns={sum(s["constant"] for s in pr["stats"])}; near-constant columns (>=95% modal value)={sum(s["near"] for s in pr["stats"])}.\n')
        f.write('  Missing columns: '+(', '.join(f'{s["name"]}={s["missing"]}' for s in pr['stats'] if s['missing']) or 'none')+'\n')
        outcols=[s for s in pr['stats'] if s['type']=='numeric']
        outtotal=sum(s['outliers'] for s in outcols)
        f.write(f'  Numeric columns={len(outcols)}; categorical columns={pr["p"]-len(outcols)}; univariate Tukey 1.5xIQR flagged cells={outtotal}. This is a screening count, not evidence of erroneous data.\n')
    # target correlations across the 1,936 numeric columns, data.csv target excluded
    y=[float(r[h2.index('vital.status')]) for r in r2]
    correlations=[]
    for col,vec in dc['vectors'].items():
        if col=='vital.status': continue
        if len(vec)==len(y): correlations.append((abs(pearson(vec,y)),pearson(vec,y),col))
    correlations.sort(reverse=True)
    f.write('\nLargest absolute Pearson correlations with numeric coding of vital.status (descriptive screening only; not causal or validated prediction evidence):\n')
    for ab,corr,col in correlations[:15]: f.write(f'  {col}: r={corr:.6f}\n')
    # screened feature intercorrelation subset: first 10 columns per group, deterministic and reproducible
    chosen=[]
    for pref in ['rs','cn','mu','pp']:
        chosen.extend([c for c in h2 if c.startswith(pref+'_')][:10])
    pairs=[]
    for i,a in enumerate(chosen):
        for b in chosen[i+1:]:
            r=pearson(dc['vectors'][a],dc['vectors'][b])
            if abs(r)>=.95: pairs.append((a,b,r))
    f.write(f'\nInter-feature correlation screen: first 10 columns in each observed prefix family (40 columns; deterministic header order); pairs with |Pearson r| >= 0.95: {len(pairs)}. This subset screen is not exhaustive over all 1,936 features.\n')
    for a,b,r in pairs[:100]: f.write(f'  {a} <> {b}: r={r:.6f}\n')
    f.write('\nSample-level QC: no exact duplicate rows in either table; data.csv has no missing cells; enriched table has missing receptor-status values but no missing feature values; no ID or sample-quality indicators were supplied. Numeric ranges and univariate outliers are summarized in dataset_structure.csv. Raw first/last row values are not reproduced; per-column SHA-256 fingerprints of first/last five entries are supplied there instead.\n')

# research readiness
with (OUT/'research_readiness_report.txt').open('w',encoding='utf-8') as f:
    f.write('RESEARCH READINESS AND METHOD SUITABILITY\n\n')
    f.write('Dataset size is 705 rows with 1,937 numeric columns in data.csv (including target), or 1,936 numeric non-target candidate features. This is p > n and has only 94 rows with target code 1. Any predictive estimate is vulnerable to high variance, overfit, and imbalance; independent external validation is essential.\n\n')
    f.write('Classical ML: possible only as exploratory, carefully nested patient-level validation after label semantics and provenance are confirmed. Logistic regression and linear SVM can serve as regularized baselines. Random forest and gradient boosting/XGBoost can be compared cautiously, but 705 samples limit reliable tuning. Feature selection must occur inside each training fold.\n')
    f.write('Multi-omics integration: feature families are co-located in one matrix, but labels/assays and patient matching are undocumented. A technical joint-matrix experiment is possible after provenance and layer semantics are established; current evidence is insufficient for a scientifically verified multimodal claim.\n')
    f.write('Deep learning: not justified for this 705-row, high-dimensional tabular dataset without additional structure or substantially larger cohort. QLattice/symbolic models may be explored after target semantics and leakage are resolved, with strong nested validation and simple baselines; small sample size and 1,935 candidate predictors create substantial selection optimism.\n')
    f.write('XAI: model equations/feature use can be inspected, but explanations are not stable or biologically interpretable until variables are documented and estimates validated.\n\n')
    f.write('Readiness: NOT ready for model development. Resolve source/provenance, row unit, target code mapping and time horizon, feature dictionary/assays, measurement timing, and patient-level IDs. Missingness is numerically low overall but concentrated in receptor-label fields. There are 99 exact duplicate-valued extra columns; investigate source provenance before any later handling. The majority/minority ratio is 6.50:1. Potential target leakage exists because vital.status is itself a predictor column and the timing of other features is unknown.\n\n')
    f.write('Source statement: the archive README links to the Kaggle dataset “Multi-Omics Breast Cancer Prediction QLattice” and describes TCGA origin only as “likely derived from TCGA or similar.” The archive does not prove the original cohort or TCGA provenance. The Kaggle page did not provide readable metadata during this audit. No provenance is inferred beyond that.\n')

with (OUT/'dataset_c_audit_complete.md').open('w',encoding='utf-8') as f:
    f.write('''# Dataset C audit and exploratory analysis\n\n## Scope and method\n\nThis report analyzes only `archive (8).zip` as Dataset C. It was extracted into `dataset_c_extracted`; no other breast cancer project or dataset was included. No model was trained, no source data were modified, and no feature/sample was removed. The archive contains two tabular datasets, a README, a notebook, and two existing QLattice model files. The audit calculations are reproducible using the included Python script and Python standard library. This environment did not have pandas, NumPy, SciPy, scikit-learn, matplotlib, seaborn, or nbformat installed; plots were generated as SVG directly.\n\n## Dimensions and files\n\n| File | Rows | Columns | Numeric | Categorical | Missing cells | Duplicate rows | Identical-value column groups |\n|---|---:|---:|---:|---:|---:|---:|---:|\n''')
    for fn,pr in profiles.items():
        ncat=sum(s['type']=='categorical' for s in pr['stats']); miss=sum(s['missing'] for s in pr['stats']); cells=pr['n']*pr['p']; dups=collections.defaultdict(list)
        for j,c in enumerate(pr['header']): dups[tuple(r[j] for r in pr['rows'])].append(c)
        dg=[x for x in dups.values() if len(x)>1]
        f.write(f'| `{fn}` | {pr["n"]} | {pr["p"]} | {pr["p"]-ncat} | {ncat} | {miss} ({100*miss/cells:.6f}%) | {pr["n"]-len({tuple(r) for r in pr["rows"]})} | {len(dg)} groups, {sum(len(x)-1 for x in dg)} extra columns |\n')
    f.write('''\nAll columns and descriptive statistics are in `dataset_structure.csv`. First/last five rows were inspected computationally; to avoid exposing sample-level records in the report, the structure CSV contains per-column SHA-256 fingerprints for each five-row edge instead of raw values.\n\n## Target and task\n\nThe only target explicitly used by notebook code and existing QLattice artifacts is `vital.status`. It is binary-coded: 0 = 611/705 (86.67%); 1 = 94/705 (13.33%); majority/minority ratio = 6.50:1. **The actual clinical meaning of 0 and 1 is not verified.** The README lists subtype or survival prediction as general possibilities, but does not define these labels. There is no survival duration, censoring, or named molecular subtype column. A provisional binary classification task is the only task evidenced by the implementation; it is not ready until the code mapping and prediction time are documented. Full counts for `PR.Status`, `ER.Status`, `HER2.Final.Status`, and `histological.type` are in `target_analysis.txt` and `CLASS` summaries in the figures.\n\n## Modalities and integration\n\nThe table has 604 `rs_`, 249 `mu_`, 860 `cn_`, and 223 `pp_` prefixed numeric columns. These are observed feature families, not confirmed modalities. Names and README context are consistent with some biological measurements, but no data dictionary establishes their assays or values. In particular, `rs_` must not be assumed to mean SNP or expression. Four categorical clinical/phenotype-like columns are available only in `brca_data_w_subtypes.csv`: PR, ER, HER2 status, and histological type. No demographic features or independent clinical metadata file are present.\n\nBoth tables have 705 rows, share 1,937 columns exactly in the same order, and match exactly row by row over those columns. The enriched file adds the four categorical fields. No explicit patient/sample identifiers exist, so patient-level matching across omics or with external cohorts cannot be measured. Feature families are technically co-stored in one numeric matrix, but shared-patient provenance is unverified.\n\n## Missingness and data quality\n\n`data.csv` has 0 missing cells. `brca_data_w_subtypes.csv` has 389 blanks out of 1,368,405 cells (0.0284%), all in `PR.Status` (122), `ER.Status` (122), and `HER2.Final.Status` (145). No exact duplicate rows were found. Each file has 62 groups of columns with identical value vectors, totaling 99 redundant extra columns. Constant and near-constant columns, all descriptive statistics, and IQR outlier screening counts are in `dataset_structure.csv` and `data_quality_report.txt`. The 40-column reproducible correlation screen is explicitly labeled non-exhaustive; target correlations cover every numeric candidate column and are descriptive only.\n\n## Provenance, target leakage, limitations\n\nThe README links the archive to a Kaggle dataset titled “Multi-Omics Breast Cancer Prediction QLattice.” It says the source is Kaggle/Computational Oncology and says the data are “likely derived from TCGA or similar.” That wording is uncertain, and no accession, cohort citation, codebook, or IDs are included. TCGA provenance is therefore not established.\n\n`vital.status` is included as a column in the candidate predictor table and must be excluded from predictors if it is the target. If receptor or histology labels become targets, related columns may leak the outcome; without temporal metadata, leakage risk cannot be resolved. Other risks include p > n, only 94 examples in the minority-coded class, 99 redundant columns, unknown sample unit, uncertain modality meanings, no external validation, and possible optimistic selection in the existing notebook. Its narrative metrics are not reproduced or independently validated here.\n\n## Readiness and next experiments\n\n**Not ready for model development.** First obtain authoritative provenance, a codebook defining `vital.status`, sample unit and prediction horizon, feature/assay dictionary, timing, and patient/sample IDs. Preserve both source tables; investigate duplicate-valued columns with the publisher. Once clarified: (1) define one clinical estimand; (2) specify patient-level nested/repeated stratified CV and an external holdout; (3) keep all imputation/scaling/selection inside training folds; (4) compare prevalence/regularized logistic and linear SVM baselines; (5) evaluate modality groups separately only after labels are verified; (6) test integration only with documented patient alignment; (7) consider tree boosting and QLattice against baselines with nested tuning; (8) report PR-AUC, ROC-AUC, sensitivity, specificity, calibration, confidence intervals, and permutation/statistical checks; (9) assess explanation stability; (10) validate externally. Deep learning is not justified by the current 705-row tabular sample. No experiments above were run.\n\n## Figures\n\nSee `figures/` for target distribution, missingness, feature family counts, example numerical distributions, a feature correlation heatmap, and row-alignment diagram. The alignment visualization depicts positional equality, not patient-ID matching.\n\n## Reproducibility\n\nRun from the workspace root: `python dataset_c_audit/run_audit.py`. It reads only `dataset_c_extracted/data.csv` and `dataset_c_extracted/brca_data_w_subtypes.csv` and rewrites audit outputs in `dataset_c_audit/`. No row-level values or identifiers are written to reports.\n''')

# Figures (SVG, deterministic)
# target distribution
body='<text x="50" y="40" class="title">Dataset C: vital.status coded classes</text>'
classes=[('0',611,'#3b82f6'),('1',94,'#f97316')]
for i,(lab,n,col) in enumerate(classes):
    x=120+i*360; h=350*n/611; y=490-h
    body+=f'<rect x="{x}" y="{y}" width="180" height="{h}" fill="{col}"/><text x="{x+70}" y="{y-12}" class="label">{n} ({100*n/705:.2f}%)</text><text x="{x+84}" y="520" class="label">{lab}</text>'
body+='<text x="380" y="570" class="small">Meaning of 0/1 is not documented in the archive</text>'
svg_write('target_distribution.svg',body)
# missingness across columns
body='<text x="40" y="36" class="title">Missing values by column (Dataset C)</text>'
all_missing=[]
for fn,pr in profiles.items():
    for s in pr['stats']:
        if s['missing']: all_missing.append((fn,s['name'],s['missing']))
body+='<text x="40" y="75" class="label">data.csv: 0 missing cells</text>'
for i,(fn,col,n) in enumerate(all_missing):
    y=120+i*70; body+=f'<text x="40" y="{y+20}" class="label">{esc(col)}</text><rect x="300" y="{y}" width="{n*3}" height="28" fill="#dc6b45"/><text x="{310+n*3}" y="{y+20}" class="label">{n} / 705</text>'
svg_write('missing_values.svg',body,height=400)
# feature categories counts
body='<text x="40" y="36" class="title">Observed feature-name families (not verified modalities)</text>'
colors=['#5b8ff9','#61d9a4','#65789b','#f6bd16']; fam=[('rs_',604),('cn_',860),('mu_',249),('pp_',223)]
for i,(lab,n) in enumerate(fam):
 y=100+i*90; body+=f'<text x="50" y="{y+25}" class="label">{lab}</text><rect x="180" y="{y}" width="{n/860*560:.1f}" height="42" fill="{colors[i]}"/><text x="{190+n/860*560:.1f}" y="{y+26}" class="label">{n}</text>'
body+='<text x="50" y="500" class="small">The archive does not define these prefixes with a feature dictionary.</text>'
svg_write('feature_family_counts.svg',body)
# numerical distributions: example first column each family, histogram 20 bins
chosen=[]
for pref in ['rs','cn','mu','pp']:
 chosen.append(next(c for c in h2 if c.startswith(pref+'_')))
body='<text x="35" y="35" class="title">Example numeric feature distributions</text>'
for panel,col in enumerate(chosen):
 vals=dc['vectors'][col]; lo=min(vals); hi=max(vals); bins=20; counts=[0]*bins
 for v in vals:
  idx=min(bins-1,int((v-lo)/(hi-lo)*bins)) if hi>lo else 0; counts[idx]+=1
 maxc=max(counts); x0=50+(panel%2)*490; y0=80+(panel//2)*240
 body+=f'<text x="{x0}" y="{y0}" class="label">{esc(col)}</text>'
 for k,n in enumerate(counts):
  bw=380/bins; hh=130*n/maxc if maxc else 0; body+=f'<rect x="{x0+k*bw}" y="{y0+170-hh}" width="{bw-2}" height="{hh}" fill="#5b8ff9"/>'
 body+=f'<text x="{x0}" y="{y0+190}" class="small">range {lo:.3g} to {hi:.3g}; n={len(vals)}</text>'
svg_write('numeric_feature_distributions.svg',body,height=560)
# correlation heatmap on first 8 of each family
selected=[]
for pref in ['rs','cn','mu','pp']: selected += [c for c in h2 if c.startswith(pref+'_')][:8]
body='<text x="30" y="28" class="title">Pearson correlation heatmap: first 8 columns per observed family</text>'
cell=15; x0=190; y0=70
for i,a in enumerate(selected):
 body+=f'<text transform="translate({x0-5},{y0+i*cell+11}) rotate(-45)" class="small">{esc(a)}</text>'
 for j,b in enumerate(selected):
  r=pearson(dc['vectors'][a],dc['vectors'][b]); shade=int(255*(1-abs(r))); color=f'rgb({255 if r<0 else shade},{shade},{255 if r>=0 else shade})'
  body+=f'<rect x="{x0+j*cell}" y="{y0+i*cell}" width="{cell}" height="{cell}" fill="{color}"/>'
svg_write('correlation_heatmap.svg',body,width=900,height=700)
# alignment diagram
body='<text x="40" y="40" class="title">Dataset C table alignment (positional evidence)</text>'
body+='<rect x="80" y="100" width="330" height="220" rx="12" fill="#e5efff" stroke="#4776b4"/><text x="120" y="135" class="label">data.csv</text><text x="110" y="175" class="label">705 rows × 1,937 columns</text><text x="110" y="215" class="small">shared numeric matrix + vital.status</text>'
body+='<rect x="590" y="100" width="330" height="220" rx="12" fill="#e6f7ed" stroke="#479365"/><text x="625" y="135" class="label">brca_data_w_subtypes.csv</text><text x="620" y="175" class="label">705 rows × 1,941 columns</text><text x="620" y="215" class="small">same 1,937 columns + 4 labels</text>'
body+='<line x1="410" y1="210" x2="590" y2="210" stroke="#333" stroke-width="3"/><polygon points="580,202 595,210 580,218" fill="#333"/><text x="433" y="190" class="small">705 / 705 rows match by position</text><text x="120" y="375" class="label">Patient ID columns found: 0</text><text x="120" y="410" class="small">Positional table alignment is not verified patient-level matching.</text>'
svg_write('sample_overlap.svg',body,height=480)

print('Audit generated:',len(list(OUT.iterdir())),'top-level entries; figures:',len(list((OUT/'figures').glob('*.svg'))))


