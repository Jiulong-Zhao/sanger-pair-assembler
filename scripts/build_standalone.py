"""Generate a self-contained, single-file CLI from the maintained package modules.

Only Biopython is required at runtime. Rebuild before releasing a new version.
"""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
src = root / 'src' / 'sanger_assembler'
version = re.search(r'__version__\s*=\s*"([^"]+)"',
                    (src / '__init__.py').read_text(encoding='utf-8')).group(1)

core = (src / 'core.py').read_text(encoding='utf-8')
cli = (src / 'cli.py').read_text(encoding='utf-8')
batch = (src / 'batch.py').read_text(encoding='utf-8')
core = core.replace('from __future__ import annotations\n', '')
cli = cli.replace('from __future__ import annotations\n', '')
batch = batch.replace('from __future__ import annotations\n', '')
cli = cli.replace('from . import __version__\n', '')
cli = cli.replace('from .core import assemble, read_sequence, reverse_complement, to_fasta\n', '')
cli = cli.replace('            from .batch import run_batch\n', '')
batch = batch.replace('from .cli import output_paths, run\n', '')
batch = batch.replace('from .core import to_fasta\n', '')
header = '''#!/usr/bin/env python3
"""Sanger Pair Assembler: standalone, batch + single-pair version.

Generated from src/sanger_assembler/{{core,cli,batch}}.py.
Use only Biopython as a dependency. MIT license: see LICENSE in the project.
"""
from __future__ import annotations
__version__ = "{version}"
'''.format(version=version)
output = header + '\n\n' + core + '\n\n' + cli + '\n\n' + batch + '''\n\nif __name__ == '__main__':
    raise SystemExit(main())
'''
result = root / 'sanger_merge.py'
result.write_text(output, encoding='utf-8')
result.chmod(0o755)
print('Generated {} ({} bytes)'.format(result, len(output.encode('utf-8'))))
