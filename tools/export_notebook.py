#!/usr/bin/env python3
"""Copy the working notebook into the repository, safely.

    python tools/export_notebook.py ../main.ipynb              # clear outputs
    python tools/export_notebook.py ../main.ipynb --keep-outputs

Executed outputs of this notebook contain participant identifiers (sub-XXXX)
and, in some cells, identifier alongside treatment arm, motion and MADRS. With
--keep-outputs those identifiers are replaced by stable pseudonyms (P01, P02,
...) unless --no-pseudonymise is given. The mapping is NOT written anywhere.

Always refuses to write a notebook that still contains an absolute path.
"""
import argparse, json, re, sys
from pathlib import Path

ID = re.compile(r'sub-(\d{3,4})')


def outputs_text(cell):
    for o in cell.get('outputs', []):
        for key in ('text',):
            if key in o:
                yield o, key, None
        for mime in ('text/plain', 'text/html'):
            if mime in o.get('data', {}):
                yield o, 'data', mime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source', type=Path, nargs='?',
                    help='working notebook; defaults to WORKING_NOTEBOOK in config_local.py')
    ap.add_argument('--keep-outputs', action='store_true')
    ap.add_argument('--no-pseudonymise', action='store_true')
    ap.add_argument('--with-results', action='store_true',
                    help='also copy SAVEFIG and SAVETABLE contents into results/')
    ap.add_argument('--dest', type=Path,
                    default=Path(__file__).resolve().parent.parent / 'notebooks' / 'main.ipynb')
    a = ap.parse_args()

    repo = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo))
    import config as cfg

    if a.source is None:
        src = getattr(cfg, 'WORKING_NOTEBOOK', None)
        if not src:
            sys.exit('no source given and WORKING_NOTEBOOK is not set in config_local.py')
        a.source = Path(src)
    if not a.source.exists():
        sys.exit(f'source notebook not found: {a.source}')

    nb = json.loads(a.source.read_text())

    abs_paths = sum(len(re.findall(r'/Users/|C:\\\\', ''.join(c['source'])))
                    for c in nb['cells'])
    if abs_paths:
        sys.exit(f'refusing: {abs_paths} absolute path(s) still in the source. '
                 'Rewire them to config.py first.')

    seen, n_ids, n_out = {}, 0, 0
    for c in nb['cells']:
        if c.get('cell_type') != 'code':
            continue
        if not a.keep_outputs:
            n_out += len(c.get('outputs', []))
            c['outputs'], c['execution_count'] = [], None
            continue
        if a.no_pseudonymise:
            continue
        for o, key, mime in outputs_text(c):
            target = o[key] if mime is None else o['data']
            k = mime if mime else None
            text = ''.join(target if k is None else target[k])

            def sub(m):
                nonlocal n_ids
                n_ids += 1
                return seen.setdefault(m.group(0), f'P{len(seen) + 1:02d}')

            new = ID.sub(sub, text)
            new = re.sub(r'/Users/[^/\s]+', '~', new)      # no usernames/paths in outputs
            new = re.sub(r'/Volumes/[^\n"]+', '<local path>', new)  # external-drive paths
            if k is None:
                o[key] = new.splitlines(keepends=True)
            else:
                o['data'][k] = new.splitlines(keepends=True)

    a.dest.write_text(json.dumps(nb, indent=1, ensure_ascii=False))

    if a.keep_outputs:
        print(f'wrote {a.dest} with outputs; '
              f'{n_ids} identifier mention(s) across {len(seen)} participants '
              f'{"left as-is" if a.no_pseudonymise else "pseudonymised"}')
        if a.no_pseudonymise:
            print('WARNING: participant identifiers are present in the committed notebook.')
    else:
        print(f'wrote {a.dest}; cleared {n_out} outputs')

    if a.with_results:
        import shutil
        for src_dir, sub in ((cfg.SAVEFIG, 'figures/' + cfg.TAG + f'_n{cfg.N_SUBJECTS}'),
                             (cfg.SAVETABLE, 'tables')):
            src_dir = Path(src_dir)
            dst = repo / 'results' / sub
            if not src_dir.is_dir():
                print(f'skipped {src_dir} (not found)'); continue
            dst.mkdir(parents=True, exist_ok=True)
            n = 0
            for f in src_dir.iterdir():
                if f.is_file() and not f.name.startswith('.') \
                        and 'per_participant' not in f.name:
                    shutil.copy2(f, dst / f.name); n += 1
            print(f'copied {n} file(s) -> results/{sub}')

    remaining = ID.findall(a.dest.read_text())
    if remaining and a.no_pseudonymise is False:
        print(f'WARNING: {len(remaining)} sub-XXXX still present after processing')


if __name__ == '__main__':
    main()
