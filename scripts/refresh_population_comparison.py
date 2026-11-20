"""Update only the mutation-specific comparison; preserve saved cohort curves."""
import csv
import gzip
import io
import json
from pathlib import Path
from build_population_coverage import draw
from refresh_dominant import save

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'dist'


def model_percentage(guides, count):
    return min(guides / count, 1) * 100 if count else None


def build():
    rendered = 0
    for path in sorted((OUT / 'data').glob('population-*.json')):
        letter = path.stem.removeprefix('population-')
        current = json.loads((OUT / 'data' / f'dominant-{letter}.json').read_text())
        shard = json.loads(path.read_text())
        for symbol, wrapper in shard.items():
            population = wrapper['population_coverage']
            dominant = current[symbol]['dominant_current']
            for strategy, series in population.get('strategies', {}).items():
                series.pop('current_comparison', None)
                if strategy == 'excision' or not dominant['assessed']:
                    continue
                name = f'{symbol}-{strategy}-current'
                plot = OUT / 'plots/population' / (name + '.svg')
                draw(symbol, strategy, {'dominant_mutations': dominant['count']}, series, population['haplotypes'], plot,
                     baseline_label=f"Mutation-specific model ({dominant['count']} mapped dominant alleles; {dominant['date']})",
                     description=f"Saved cohort curve; current ClinVar mapped subset {dominant['count']} alleles as of {dominant['date']}; {dominant['unresolved']} P/LP alleles with unresolved inheritance excluded. Equal-frequency model, not observed disease prevalence.")
                buffer = io.StringIO()
                writer = csv.writer(buffer)
                writer.writerow(['prioritized_gRNAs', 'new_haplotypes', 'cumulative_haplotypes', 'cohort_haplotypes', 'coverage_percent', 'prePAM_possible_haplotypes', 'selected_allele_1', 'selected_allele_2', 'mutation_specific_model_percent', 'mapped_dominant_alleles', 'inheritance_unresolved_alleles', 'ClinVar_release'])
                for point in [dict(guides=0, added=0, covered=0, allele1='', allele2='')] + series['points']:
                    writer.writerow([point['guides'], point['added'], point['covered'], population['haplotypes'], point['covered']/population['haplotypes']*100, series['possible'], point['allele1'], point['allele2'], model_percentage(point['guides'], dominant['count']), dominant['count'], dominant['unresolved'], dominant['date']])
                export = OUT / 'downloads/population' / (name + '.csv.gz')
                export.write_bytes(gzip.compress(buffer.getvalue().encode(), mtime=0))
                series['current_comparison'] = dict(count=dominant['count'], unresolved=dominant['unresolved'], date=dominant['date'], plot=str(plot.relative_to(OUT)), export=dict(path=str(export.relative_to(OUT))))
                rendered += 1
                if rendered % 100 == 0:
                    print(f'Updated {rendered} population comparisons', flush=True)
        save(path, shard)
    print(f'Updated {rendered} comparisons; saved cohort coverage and original plots preserved.', flush=True)


if __name__ == '__main__':
    build()
