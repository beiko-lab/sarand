import json
import re
from pathlib import Path


def read_alignment_deltas(root_dir):
    out = dict()
    align_dir = root_dir / 'AMR_info' / 'alignments'
    for file in align_dir.glob('*_delta.json'):
        with open(file, 'r') as f:
            data = json.load(f)
            key = re.match(r'(amr_group_\d+).+', file.name).group(1)
            assert(key not in out)
            out[key] = data
    return out


def report_on_align_deltas(deltas, root_dir):
    header = ['amr_group', 'amr', 'result', 'start_pos', 'end_pos', 'nodes', 'orientations']
    lines = ['\t'.join(header)]
    for amr_group, d_delta in sorted(deltas.items(), key=lambda x: int(x[0].split('_')[-1])):
        ba = d_delta['bandage']
        ga = d_delta['graph_aligner']

        keys_in_ba = set(ba.keys())
        keys_in_ga = set(ga.keys())

        keys_in_both = keys_in_ba.intersection(keys_in_ga)
        keys_in_ba_only = keys_in_ba.difference(keys_in_ga)
        keys_in_ga_only = keys_in_ga.difference(keys_in_ba)

        for amr in sorted(keys_in_ba_only):
            for d_amr in ba[amr]:
                row = [
                    amr_group,
                    amr,
                    'amr_in_bandage_only',
                    d_amr['start_pos'],
                    d_amr['end_pos'],
                    d_amr['nodes'],
                    d_amr['orientations']
                ]
                lines.append('\t'.join(map(str, row)))

        for amr in sorted(keys_in_ga_only):
            for d_amr in ga[amr]:
                row = [
                    amr_group,
                    amr,
                    'amr_in_ga_only',
                    d_amr['start_pos'],
                    d_amr['end_pos'],
                    d_amr['nodes'],
                    d_amr['orientations']
                ]
                lines.append('\t'.join(map(str, row)))

        for amr in sorted(keys_in_both):
            lst_amr_ba = ba[amr]
            lst_amr_ga = ga[amr]

            # Convert the rows into strings for comparison
            set_amr_ba = convert_lst_amr_to_set(lst_amr_ba)
            set_amr_ga = convert_lst_amr_to_set(lst_amr_ga)

            set_items_in_ba_only = set_amr_ba.difference(set_amr_ga)
            set_items_in_ga_only = set_amr_ga.difference(set_amr_ba)
            set_items_in_common = set_amr_ba.intersection(set_amr_ga)

            for item in set_items_in_ba_only:
                row = [
                    amr_group,
                    amr,
                    'bandage_different',
                    item[0],
                    item[1],
                    item[2],
                    item[3]
                ]
                lines.append('\t'.join(map(str, row)))

            for item in set_items_in_ga_only:
                row = [
                    amr_group,
                    amr,
                    'ga_different',
                    item[0],
                    item[1],
                    item[2],
                    item[3]
                ]
                lines.append('\t'.join(map(str, row)))

            for item in set_items_in_common:
                row = [
                    amr_group,
                    amr,
                    'same',
                    item[0],
                    item[1],
                    item[2],
                    item[3]
                ]
                lines.append('\t'.join(map(str, row)))

    with open(root_dir / 'align_deltas.tsv', 'w') as f:
        f.write('\n'.join(lines))
    return


def convert_lst_amr_to_set(lst):
    out = set()
    for item in lst:
        cur_item = [
            item['start_pos'],
            item['end_pos'],
            str(item['nodes']),
            str(item['orientations'])
        ]
        out.add(tuple(cur_item))
    return out


def read_annotation_deltas(root_dir):
    out = dict()
    anno_dir = root_dir / 'annotations' / 'annotations_1000'
    for delta_path in anno_dir.glob('*/delta_*.json'):
        with open(delta_path, 'r') as f:
            data = json.load(f)
            amr_name = re.match(r'annotation_(.+)_1000', delta_path.parent.name).group(1)
            amr_name = f"{amr_name}_{delta_path.stem.split('_')[1]}"
            assert(amr_name not in out)
            for x in data['prokka']:
                del x['prokka_gene_name']
            for x in data['bakta']:
                del x['prokka_gene_name']
            out[amr_name] = data
    return out


def convert_anno_to_keyed(d_anno):
    out = dict()
    for item in d_anno:
        key = (item['start_pos'], item['end_pos'])
        assert (key not in out)
        out[key] = item
    return out


def report_on_anno_deltas(deltas, root_dir):
    header = [
        'amr', 'uid', 'source', 'result', 'locus_tag', 'gene', 'length', 'product', 'start_pos',
        'end_pos', 'RGI_prediction_type', 'coverage',
        'family', 'seq_name', 'target_amr', 'seq_value'
    ]
    lines = ['\t'.join(header)]

    result_id = 0

    for amr_name, d_delta in sorted(deltas.items()):
        d_prokka = d_delta['prokka']
        d_bakta = d_delta['bakta']

        d_prokka_keyed = convert_anno_to_keyed(d_prokka)
        d_bakta_keyed = convert_anno_to_keyed(d_bakta)

        prokka_keys = set(d_prokka_keyed.keys())
        bakta_keys = set(d_bakta_keyed.keys())

        keys_in_both = prokka_keys.intersection(bakta_keys)
        keys_in_prokka_only = prokka_keys.difference(bakta_keys)
        keys_in_bakta_only = bakta_keys.difference(prokka_keys)

        for key in sorted(keys_in_prokka_only):
            cur_dict = d_prokka_keyed[key]
            cur_row = [
                amr_name,
                '',
                'prokka',
                'not_in_bakta',
                cur_dict['locus_tag'],
                cur_dict['gene'],
                cur_dict['length'],
                cur_dict['product'],
                cur_dict['start_pos'],
                cur_dict['end_pos'],
                cur_dict['RGI_prediction_type'],
                cur_dict['coverage'],
                cur_dict['family'],
                cur_dict['seq_name'],
                cur_dict['target_amr'],
                cur_dict['seq_value']
            ]
            lines.append('\t'.join(map(str, cur_row)))

        for key in sorted(keys_in_bakta_only):
            cur_dict = d_bakta_keyed[key]
            cur_row = [
                amr_name,
                '',
                'bakta',
                'not_in_prokka',
                cur_dict['locus_tag'],
                cur_dict['gene'],
                cur_dict['length'],
                cur_dict['product'],
                cur_dict['start_pos'],
                cur_dict['end_pos'],
                cur_dict['RGI_prediction_type'],
                cur_dict['coverage'],
                cur_dict['family'],
                cur_dict['seq_name'],
                cur_dict['target_amr'],
                cur_dict['seq_value']
            ]
            lines.append('\t'.join(map(str, cur_row)))

        for key in sorted(keys_in_both):
            cur_prokka = d_prokka_keyed[key]
            cur_bakta = d_bakta_keyed[key]

            if cur_prokka == cur_bakta:
                cur_row = [
                    amr_name,
                    '',
                    'both',
                    'same',
                    cur_bakta['locus_tag'],
                    cur_bakta['gene'],
                    cur_bakta['length'],
                    cur_bakta['product'],
                    cur_bakta['start_pos'],
                    cur_bakta['end_pos'],
                    cur_bakta['RGI_prediction_type'],
                    cur_bakta['coverage'],
                    cur_bakta['family'],
                    cur_bakta['seq_name'],
                    cur_bakta['target_amr'],
                    cur_bakta['seq_value']
                ]
                lines.append('\t'.join(map(str, cur_row)))

            else:
                cur_row = [
                    amr_name,
                    result_id,
                    'prokka',
                    'different',
                    cur_prokka['locus_tag'],
                    cur_prokka['gene'],
                    cur_prokka['length'],
                    cur_prokka['product'],
                    cur_prokka['start_pos'],
                    cur_prokka['end_pos'],
                    cur_prokka['RGI_prediction_type'],
                    cur_prokka['coverage'],
                    cur_prokka['family'],
                    cur_prokka['seq_name'],
                    cur_prokka['target_amr'],
                    cur_prokka['seq_value']
                ]
                lines.append('\t'.join(map(str, cur_row)))

                cur_row = [
                    amr_name,
                    result_id,
                    'bakta',
                    'different',
                    cur_bakta['locus_tag'],
                    cur_bakta['gene'],
                    cur_bakta['length'],
                    cur_bakta['product'],
                    cur_bakta['start_pos'],
                    cur_bakta['end_pos'],
                    cur_bakta['RGI_prediction_type'],
                    cur_bakta['coverage'],
                    cur_bakta['family'],
                    cur_bakta['seq_name'],
                    cur_bakta['target_amr'],
                    cur_bakta['seq_value']
                ]
                lines.append('\t'.join(map(str, cur_row)))

                cur_row = [
                    amr_name,
                    result_id,
                    'delta',
                    '',
                    '' if cur_prokka['locus_tag'] == cur_bakta['locus_tag'] else '^^^^^^',
                    '' if cur_prokka['gene'] == cur_bakta['gene'] else '^^^^^^',
                    '' if cur_prokka['length'] == cur_bakta['length'] else '^^^^^^',
                    '' if cur_prokka['product'] == cur_bakta['product'] else '^^^^^^',
                    '' if cur_prokka['start_pos'] == cur_bakta['start_pos'] else '^^^^^^',
                    '' if cur_prokka['end_pos'] == cur_bakta['end_pos'] else '^^^^^^',
                    '' if cur_prokka['RGI_prediction_type'] == cur_bakta['RGI_prediction_type'] else '^^^^^^',
                    '' if cur_prokka['coverage'] == cur_bakta['coverage'] else '^^^^^^',
                    '' if cur_prokka['family'] == cur_bakta['family'] else '^^^^^^',
                    '' if cur_prokka['seq_name'] == cur_bakta['seq_name'] else '^^^^^^',
                    '' if cur_prokka['target_amr'] == cur_bakta['target_amr'] else '^^^^^^',
                    '' if cur_prokka['seq_value'] == cur_bakta['seq_value'] else '^^^^^^',
                ]
                lines.append('\t'.join(map(str, cur_row)))

                result_id += 1

    with open(root_dir / 'anno_deltas.tsv', 'w') as f:
        f.write('\n'.join(lines))

    return


def read_identity_cov_deltas(root_dir):
    out = {
        'bandage': dict(),
        'graph_aligner': dict()
    }
    align_dir = root_dir / 'AMR_info' / 'alignments'
    for file in align_dir.glob('*_coverage_identity.json'):
        with open(file, 'r') as f:
            data = json.load(f)
            if 'bandage' in file.name:
                program = 'bandage'
            elif 'graphaligner' in file.name:
                program = 'graph_aligner'
            else:
                raise Exception('Unknown program')
            for row in data:
                path_hit = re.match(r'\((\d+)\)?.+\((\d+)\)?', row['path'])
                if path_hit is None:
                    path_start = row['path']
                    path_end = row['path']
                else:
                    path_start = int(path_hit.group(1))
                    path_end = int(path_hit.group(2))

                    if program == 'bandage':
                        path_start -= 1  # due to indexing

                row_id = f'{row["id"]}_{path_start}_{path_end}'
                if row_id in out[program]:
                    row_id += '_2'
                assert(row_id not in out[program])
                out[program][row_id] = row
    return out

def key_to_path_start_end(key):
    hit = re.match(r'.+_(\d+)_(\d+)$', key)
    if not hit:
        return 'N/A', 'N/A'
    return int(hit.group(1)), int(hit.group(2))

def is_discard(iden, cov):
    if iden >= 95 and cov >= 95:
        return 'false'
    else:
        return 'true'

def report_on_identity_cov_deltas(deltas, out_dir):
    out_path = out_dir / 'identity_cov_deltas.tsv'
    header = ['id', 'discard', 'from_program', 'status', 'path_start', 'path_end', 'identity', 'coverage']
    lines = ['\t'.join(header)]

    bandage_keys = set(deltas['bandage'].keys())
    graph_aligner_keys = set(deltas['graph_aligner'].keys())

    keys_not_in_bandage = graph_aligner_keys - bandage_keys
    keys_not_in_graph_aligner = bandage_keys - graph_aligner_keys
    keys_in_both = bandage_keys & graph_aligner_keys

    for key in sorted(keys_not_in_bandage):
        path_start, path_end = key_to_path_start_end(key)
        cur_row = [
            key,
            is_discard(deltas['graph_aligner'][key]['identity'], deltas['graph_aligner'][key]['coverage']),
            'graph_aligner',
            'not_in_bandage',
            path_start,
            path_end,
            deltas['graph_aligner'][key]['identity'],
            deltas['graph_aligner'][key]['coverage']
        ]
        lines.append('\t'.join(map(str, cur_row)))

    for key in sorted(keys_not_in_graph_aligner):
        path_start, path_end = key_to_path_start_end(key)
        cur_row = [
            key,
            is_discard(deltas['bandage'][key]['identity'], deltas['bandage'][key]['coverage']),
            'bandage',
            'not_in_graph_aligner',
            path_start,
            path_end,
            deltas['bandage'][key]['identity'],
            deltas['bandage'][key]['coverage']
        ]
        lines.append('\t'.join(map(str, cur_row)))

    for key in sorted(keys_in_both):

        path_start, path_end = key_to_path_start_end(key)

        bdg_identity = round(deltas['bandage'][key]['identity'], 3)
        ga_identity = round(deltas['graph_aligner'][key]['identity'], 3)

        bdg_cov = round(deltas['bandage'][key]['coverage'], 3)
        ga_cov = round(deltas['graph_aligner'][key]['coverage'], 3)

        if bdg_identity == ga_identity and bdg_cov == ga_cov:
            cur_row = [
                key,
                is_discard(bdg_identity, bdg_cov),
                'both',
                'same',
                path_start,
                path_end,
                bdg_identity,
                bdg_cov
            ]
            lines.append('\t'.join(map(str, cur_row)))

        else:
            cur_row = [
                key,
                is_discard(bdg_identity, bdg_cov),
                'bandage',
                'different',
                path_start,
                path_end,
                bdg_identity,
                bdg_cov
            ]
            lines.append('\t'.join(map(str, cur_row)))

            cur_row = [
                key,
                is_discard(ga_identity, ga_cov),
                'graph_aligner',
                'different',
                path_start,
                path_end,
                ga_identity,
                ga_cov
            ]
            lines.append('\t'.join(map(str, cur_row)))

            cur_row = [
                key,
                '',
                'delta',
                '',
                '',
                '',
                '' if bdg_identity == ga_identity else '^^^^^^',
                '' if bdg_cov == ga_cov else '^^^^^^'
            ]
            lines.append('\t'.join(map(str, cur_row)))

    print(out_path)
    with out_path.open('w') as f:
        f.write('\n'.join(lines))
    return


def main():
    root_dir = Path('/Users/aaron/git/sarand/testing/new/1_ERR1713331/output')
    # root_dir = Path('/private/tmp/sarand_out/out')

    iden_cov_deltas = read_identity_cov_deltas(root_dir)
    report_on_identity_cov_deltas(iden_cov_deltas, root_dir)

    aln_deltas = read_alignment_deltas(root_dir)
    report_on_align_deltas(aln_deltas, root_dir)

    anno_deltas = read_annotation_deltas(root_dir)
    report_on_anno_deltas(anno_deltas, root_dir)


if __name__ == '__main__':
    main()
