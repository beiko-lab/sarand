import json
import os
import re
from collections import defaultdict

from Bio import SeqIO


class OutputFile:

    def __init__(self, path):
        self.path = path
        self.data = defaultdict(list)

    def key_same(self, file, key):
        self.data[file].append(('SAME', key, '', ''))

    def key_different(self, file, key, prev, new):
        self.data[file].append(('DIFFERENT', key, prev, new))

    def key_missing(self, file, key, program):
        self.data[file].append((f'missing_from_{program}', key, '', ''))

    def write(self):
        with open(self.path, 'w') as f:
            f.write('\t'.join(['file', 'type', 'key', 'previous', 'new']) + '\n')
            for file, values in sorted(self.data.items()):
                for value in values:
                    to_write = '\t'.join(list(map(str, value)))
                    f.write(f'{file}\t{to_write}\n')


def path_to_dict(path):
    d = {'name': os.path.basename(path)}
    if os.path.isdir(path):
        d['type'] = "directory"
        d['children'] = [path_to_dict(os.path.join(path, x)) for x in os.listdir(path)]
    else:
        d['type'] = "file"
    return d


def read_amr_info_overlaps(root_dir):
    out = list()
    with open(os.path.join(root_dir, 'AMR_info', 'overlaps.txt')) as f:
        for line in f.readlines():
            hit = line.strip()
            out.append(hit)
    return out


def delim_to_json(path, sep='\t'):
    out = list()
    with open(path) as f:
        headers = [x.strip() for x in f.readline().split(sep)]
        for line in f.readlines():
            cols = line.strip().split(sep)
            out.append(dict(zip(headers, cols)))
    return out


def read_amr_info_alignments(root_dir):
    out = dict()

    re_hit = re.compile(r'(amr_group_\d+)_.*\.tsv')
    for file in os.listdir(os.path.join(root_dir, 'AMR_info', 'alignments')):
        cur_hit = re_hit.match(file)
        if cur_hit:
            path = os.path.join(root_dir, 'AMR_info', 'alignments', file)
            name = cur_hit.group(1)
            out[name] = delim_to_json(path)
    return out


def read_amr_info_sequences(root_dir):
    out = dict()
    for file in os.listdir(os.path.join(root_dir, 'AMR_info', 'sequences')):
        path = os.path.join(root_dir, 'AMR_info', 'sequences', file)
        x = SeqIO.parse(path, 'fasta')
        out[file.replace('.fasta', '')] = {v.description: str(v.seq) for k, v in SeqIO.to_dict(x).items()}
    return out


def compare_amr_info(root_dir):
    out = dict()
    # out['alignments'] = read_amr_info_alignments(root_dir)
    out['overlaps.txt'] = read_amr_info_overlaps(root_dir)
    out['sequences'] = read_amr_info_sequences(root_dir)
    return out


def parse_anno_dir(path):
    out = dict()

    cur_amr_re = re.compile(r'annotation_(.+)_1000')
    cur_amr_name = cur_amr_re.match(os.path.basename(path)).group(1)

    anno_detail_path = os.path.join(path, f'annotation_detail_{cur_amr_name}.csv')
    coverage_anno = os.path.join(path, f'coverage_annotation_30_{cur_amr_name}.csv')
    seq_comparison_path = os.path.join(path, f'seq_comparison_genes_{cur_amr_name}.txt')
    trimmed_annno_path = os.path.join(path, f'trimmed_annotation_info_{cur_amr_name}.csv')

    out['annotation_detail'] = delim_to_json(anno_detail_path, ',')
    out['coverage_annotation'] = delim_to_json(coverage_anno, ',')
    out['seq_comparison_genes'] = list()
    with open(seq_comparison_path) as f:
        for line in f.readlines():
            out['seq_comparison_genes'].append(line.strip())
    out['trimmed_annotation_info'] = delim_to_json(trimmed_annno_path, ',')

    rgi_dir = os.path.join(path, 'rgi_dir')
    rgi_json = [os.path.join(rgi_dir, x) for x in os.listdir(rgi_dir) if x.endswith('.json')][0]
    with open(rgi_json) as f:
        out['rgi_json'] = json.load(f)
    return out


def compare_annotations(root_dir):
    # Is there any other possibility?
    files_in_anno_dir = set()
    for file in os.listdir(os.path.join(root_dir, 'annotations')):
        if os.path.isdir(os.path.join(root_dir, 'annotations', file)):
            files_in_anno_dir.add(file)
    assert (files_in_anno_dir == {'annotations_1000'})

    out = dict()

    # Read the not found in graph file
    anno_1000 = os.path.join(root_dir, 'annotations', 'annotations_1000')
    path_not_found = os.path.join(anno_1000, 'not_found_annotation_amrs_in_graph.txt')
    with open(path_not_found) as f:
        out['not_found_annotation_amrs_in_graph'] = [x.strip() for x in f.readlines()]

    # Parse each directory
    for cur_d in os.listdir(anno_1000):
        cur_d_path = os.path.join(anno_1000, cur_d)
        if os.path.isdir(cur_d_path) and not cur_d_path.endswith('_old'):
            cur_d_name = cur_d.replace('annotation_', '')
            out[cur_d_name] = parse_anno_dir(cur_d_path)
    return out


def compare_sequences(root_dir):
    # Is there any other possibility?
    files_in_anno_dir = set()
    for file in os.listdir(os.path.join(root_dir, 'sequences_info')):
        if os.path.isdir(os.path.join(root_dir, 'sequences_info', file)):
            files_in_anno_dir.add(file)
    assert (files_in_anno_dir == {'sequences_info_1000'})

    out = {
        'paths_info': dict(),
        'sequences': dict()
    }

    # paths_info
    re_paths_info_file = re.compile(r'ng_sequences_(.+)_1000_\d{4}.+\.csv')
    for file in os.listdir(os.path.join(root_dir, 'sequences_info', 'sequences_info_1000', 'paths_info')):
        path = os.path.join(root_dir, 'sequences_info', 'sequences_info_1000', 'paths_info', file)
        re_hit = re_paths_info_file.match(file)
        if re_hit:
            name = re_hit.group(1)
            out['paths_info'][name] = delim_to_json(path, ',')

    # sequences
    re_sequences_file = re.compile(r'ng_sequences_(.+)_1000_\d{4}.+\.txt')
    re_seq_file = re.compile(r'(.+):\n(.+)')
    for file in os.listdir(os.path.join(root_dir, 'sequences_info', 'sequences_info_1000', 'sequences')):
        path = os.path.join(root_dir, 'sequences_info', 'sequences_info_1000', 'sequences', file)
        re_hit = re_sequences_file.match(file)
        if re_hit:
            name = re_hit.group(1)
            with open(path) as f:
                hits = re_seq_file.findall(f.read())
                out['sequences'][name] = hits
    return out


def new_sec(title):
    return [
        '',
        '>' * 80,
        f'  {title}  '.center(80),
        '<' * 80,
    ]


def delta_amr_info_overlaps(previous_json, new_json, out_file: OutputFile):
    out = list()
    file_name = 'AMR_info/overlaps.txt'
    out.extend(new_sec(file_name))
    if new_json['amr_info']['overlaps.txt'] == previous_json['amr_info']['overlaps.txt']:
        out.append('SAME'.center(80))
        out_file.key_same(file_name, '')
    else:
        out.append(' DIFFERENT '.center(80, '@'))
        out.append(f'Previous: {previous_json["amr_info"]["overlaps.txt"]}')
        out.append(f'New:      {new_json["amr_info"]["overlaps.txt"]}')
        out_file.key_different(file_name, '', previous_json['amr_info']['overlaps.txt'],
                               new_json['amr_info']['overlaps.txt'])
    out.append('-' * 80)
    return out


def delta_anno_seqs(previous_json, new_json, out_file: OutputFile):
    prev_keys = set(previous_json['amr_info']['sequences'].keys())
    new_keys = set(new_json['amr_info']['sequences'].keys())

    prev_keys_not_in_new = prev_keys - new_keys
    new_keys_not_in_prev = new_keys - prev_keys
    common_keys = prev_keys & new_keys

    file_name = 'AMR_info/sequences'

    for missing_key in prev_keys_not_in_new:
        out_file.key_missing(file_name, f'AMR_info/sequences/{missing_key}.fasta', 'new')
    for missing_key in new_keys_not_in_prev:
        out_file.key_missing(file_name, f'AMR_info/sequences/{missing_key}.fasta', 'previous')

    for common_key in common_keys:
        prev_dict = previous_json['amr_info']['sequences'][common_key]
        new_dict = new_json['amr_info']['sequences'][common_key]
        out_str = f'AMR_info/sequences/{common_key}.fasta'
        if prev_dict == new_dict:
            out_file.key_same(out_str, '')
        else:
            out_file.key_different(out_str, '', prev_dict, new_dict)
    return list()


def delta_amr_info_alignments(previous_json, new_json, out_file: OutputFile):
    prev_keys = set(previous_json['amr_info']['alignments'].keys())
    new_keys = set(new_json['amr_info']['alignments'].keys())

    prev_keys_not_in_new = prev_keys - new_keys
    new_keys_not_in_prev = new_keys - prev_keys
    common_keys = prev_keys & new_keys

    for missing_key in prev_keys_not_in_new:
        out_file.key_missing(f'AMR_info/alignments/{missing_key}', 'FILE', 'new')
    for missing_key in new_keys_not_in_prev:
        out_file.key_missing(f'AMR_info/alignments/{missing_key}', 'FILE', 'previous')
    for common_key in common_keys:
        prev_dict = previous_json['amr_info']['alignments'][common_key]
        new_dict = new_json['amr_info']['alignments'][common_key]
        out_str = f'AMR_info/alignments/{common_key}.fasta'
        if prev_dict != new_dict:
            out_file.key_different(out_str, common_key, prev_dict, new_dict)
        else:
            out_file.key_same(out_str, '')

    return list()


def delta_anno_not_found(previous_json, new_json, out_file: OutputFile):
    file_name = 'annotations/not_found_annotation_amrs_in_graph.txt'
    if new_json['annotations']['not_found_annotation_amrs_in_graph'] == previous_json['annotations'][
        'not_found_annotation_amrs_in_graph']:
        out_file.key_same(file_name, '')
    else:
        out_file.key_different(file_name, '', previous_json['annotations']['not_found_annotation_amrs_in_graph'],
                               new_json['annotations']['not_found_annotation_amrs_in_graph'])
    return list()


def delta_anno_amr(previous_json, new_json, out_file: OutputFile):
    prev_keys = set(previous_json['annotations'].keys()) - {'not_found_annotation_amrs_in_graph'}
    new_keys = set(new_json['annotations'].keys()) - {'not_found_annotation_amrs_in_graph'}

    prev_keys_not_in_new = prev_keys - new_keys
    new_keys_not_in_prev = new_keys - prev_keys
    common_keys = prev_keys & new_keys

    for missing_key in prev_keys_not_in_new:
        out_file.key_missing(f'annotations/{missing_key}', 'FILE', 'new')
    for missing_key in new_keys_not_in_prev:
        out_file.key_missing(f'annotations/{missing_key}', 'FILE', 'previous')

    for common_key in common_keys:
        prev_dict = previous_json['annotations'][common_key]
        new_dict = new_json['annotations'][common_key]

        out_str = f'annotations/{common_key}/annotation_detail.csv'
        prev_anno_detail = prev_dict['annotation_detail']
        new_anno_detail = new_dict['annotation_detail']
        if prev_anno_detail == new_anno_detail:
            out_file.key_same(out_str, '')
        else:
            if len(new_anno_detail) == len(prev_anno_detail):
                for cur_prev_anno_detail, cur_new_anno_detail in zip(prev_anno_detail, new_anno_detail):
                    cur_prev_anno_keys = cur_prev_anno_detail.keys()
                    for cur_now_anno_key in cur_prev_anno_keys:
                        cur_prev_anno_val = cur_prev_anno_detail[cur_now_anno_key]
                        cur_new_anno_val = cur_new_anno_detail[cur_now_anno_key]
                        if cur_prev_anno_val != cur_new_anno_val:
                            out_file.key_different(out_str, cur_now_anno_key, cur_prev_anno_val, cur_new_anno_val)
            else:
                print('TODO')
                out_file.key_different(out_str, 'raw', str(prev_anno_detail), str(new_anno_detail))

        out_str = f'annotations/{common_key}/coverage_annotation.csv'
        prev_info = prev_dict['coverage_annotation']
        new_info = new_dict['coverage_annotation']
        if prev_info == new_info:
            out_file.key_same(out_str, '')
        else:
            if len(prev_info) == len(new_info):
                for cur_prev_key, cur_new_key in zip(prev_info, new_info):
                    cur_prev_keys = cur_prev_key.keys()
                    for cur_now_key in cur_prev_keys:
                        cur_prev_val = cur_prev_key[cur_now_key]
                        cur_new_val = cur_new_key[cur_now_key]
                        if cur_prev_val != cur_new_val:
                            out_file.key_different(out_str, cur_now_key, cur_prev_val, cur_new_val)
            else:
                print('TODO')
                out_file.key_different(out_str, 'raw', str(prev_info), str(new_info))

        out_str = f'annotations/{common_key}/seq_comparison_genes.txt'
        prev_info = prev_dict['seq_comparison_genes']
        new_info = new_dict['seq_comparison_genes']
        if prev_info == new_info:
            out_file.key_same(out_str, '')
        else:
            out_file.key_different(out_str, '', str(prev_info), str(new_info))

        out_str = f'annotations/{common_key}/trimmed_annotation_info.csv'
        prev_info = prev_dict['trimmed_annotation_info']
        new_info = new_dict['trimmed_annotation_info']
        if prev_info == new_info:
            out_file.key_same(out_str, '')
        else:
            if len(prev_info) == len(new_info):
                for cur_prev_key, cur_new_key in zip(prev_info, new_info):
                    cur_prev_keys = cur_prev_key.keys()
                    for cur_now_key in cur_prev_keys:
                        cur_prev_val = cur_prev_key[cur_now_key]
                        cur_new_val = cur_new_key[cur_now_key]
                        if cur_prev_val != cur_new_val:
                            out_file.key_different(out_str, cur_now_key, cur_prev_val, cur_new_val)
            else:
                print('TODO')
                out_file.key_different(out_str, 'raw', str(prev_info), str(new_info))

    return list()


def delta_sequences(previous_json, new_json, out_file: OutputFile):
    assert ({'paths_info', 'sequences'} == set(previous_json['sequences'].keys()) == set(new_json['sequences'].keys()))

    # Paths info
    prev_paths_info = previous_json['sequences']['paths_info']
    new_paths_info = new_json['sequences']['paths_info']

    prev_keys = set(prev_paths_info.keys())
    new_keys = set(new_paths_info.keys())

    prev_keys_not_in_new = prev_keys - new_keys
    new_keys_not_in_prev = new_keys - prev_keys
    common_keys = prev_keys & new_keys

    for missing_key in prev_keys_not_in_new:
        out_file.key_missing(f'sequences/paths_info/{missing_key}', 'FILE', 'new')
    for missing_key in new_keys_not_in_prev:
        out_file.key_missing(f'sequences/paths_info/{missing_key}', 'FILE', 'previous')

    for common_key in common_keys:
        prev_dict = prev_paths_info[common_key]
        new_dict = new_paths_info[common_key]

        if prev_dict == new_dict:
            out_file.key_same(f'sequences/paths_info/{common_key}', '')

        if len(prev_dict) == len(new_dict):
            for cur_prev_dict, cur_new_dict in zip(prev_dict, new_dict):
                keys = set(cur_prev_dict.keys()) & set(cur_new_dict.keys())
                for key in keys:
                    cur_prev_value = cur_prev_dict[key]
                    cur_new_value = cur_new_dict[key]
                    if cur_prev_value != cur_new_value:
                        out_file.key_different(f'sequences/paths_info/{common_key}', key, cur_prev_value, cur_new_value)
        else:
            print('TODO')
            out_file.key_different(f'sequences/paths_info/{common_key}', 'raw', str(prev_dict), str(new_dict))

    # Sequences

    prev_sequences = previous_json['sequences']['sequences']
    new_sequences = new_json['sequences']['sequences']

    prev_keys = set(prev_sequences.keys())
    new_keys = set(new_sequences.keys())

    prev_keys_not_in_new = prev_keys - new_keys
    new_keys_not_in_prev = new_keys - prev_keys
    common_keys = prev_keys & new_keys

    for missing_key in prev_keys_not_in_new:
        out_file.key_missing(f'sequences/sequences/{missing_key}', 'FILE', 'new')
    for missing_key in new_keys_not_in_prev:
        out_file.key_missing(f'sequences/sequences/{missing_key}', 'FILE', 'previous')

    for common_key in common_keys:
        prev_dict = prev_sequences[common_key]
        new_dict = new_sequences[common_key]

        if prev_dict == new_dict:
            out_file.key_same(f'sequences/sequences/{common_key}', '')
        else:
            out_file.key_different(f'sequences/sequences/{common_key}', '', str(prev_dict), str(new_dict))
    return


def do_comparison(previous_json, new_json, out_file):
    # AMR_info_alignments
    # delta_amr_info_alignments(previous_json, new_json, out_file)

    # AMR_info/overlaps.txt
    delta_amr_info_overlaps(previous_json, new_json, out_file)

    # AMR_info/sequences
    delta_anno_seqs(previous_json, new_json, out_file)

    delta_anno_not_found(previous_json, new_json, out_file)

    delta_anno_amr(previous_json, new_json, out_file)

    delta_sequences(previous_json, new_json, out_file)

    return list()


def main():
    previous_dir = '/Users/aaron/git/sarand/test/expected_output'
    new_dir = '/private/tmp/sarand_out/out'
    new_dir = '/private/tmp/sarand/1_1_1/sarand_results_docker'

    previous_dir = '/Users/aaron/git/sarand/testing/expected/CAMI_M_2/sarand_results'
    new_dir = '/Users/aaron/git/sarand/testing/new/CAMI_M_2/output'

    previous_json = {
        'amr_info': compare_amr_info(previous_dir),
        'annotations': compare_annotations(previous_dir),
        'sequences': compare_sequences(previous_dir)
    }
    new_json = {
        'amr_info': compare_amr_info(new_dir),
        'annotations': compare_annotations(new_dir),
        'sequences': compare_sequences(new_dir)
    }

    dir_out = '/tmp'
    os.makedirs(dir_out, exist_ok=True)

    # Write the raw data
    with open(os.path.join(dir_out, 'raw_data.json'), 'w') as f:
        json.dump({'previous': previous_json, 'new': new_json}, f, indent=4)

    with open(os.path.join(dir_out, 'previous.json'), 'w') as f:
        json.dump(previous_json, f, indent=4)
    with open(os.path.join(dir_out, 'new.json'), 'w') as f:
        json.dump(new_json, f, indent=4)

    delta_path = os.path.join(dir_out, 'delta.tsv')
    out_file = OutputFile(delta_path)
    do_comparison(previous_json, new_json, out_file)
    out_file.write()

    return


if __name__ == '__main__':
    main()
