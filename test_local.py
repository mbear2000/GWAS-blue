"""Offline tests; never opens SecureCRT, connects SSH, or submits real jobs."""
import base64
import os
import re
from pathlib import Path
import subprocess
import tempfile
import unittest
import server

BASH = r'C:\Program Files\Git\bin\bash.exe'


class ValidationTests(unittest.TestCase):
    def test_unique_local_logs(self):
        with tempfile.TemporaryDirectory(dir=server.ROOT) as d:
            original=Path(d)/'LSH2022_Re-fitting_GWASrun.log'
            self.assertEqual(server.available_log_path(original),original)
            original.write_text('old log','utf-8')
            second=server.available_log_path(original)
            second.write_text('second log','utf-8')
            third=server.available_log_path(original)
            self.assertEqual(second.parent,original.parent)
            self.assertNotEqual(second,third)
            self.assertFalse(third.exists())
            self.assertEqual(original.read_text(),'old log')

    def test_incremental_log_reader(self):
        with tempfile.TemporaryDirectory(dir=server.ROOT) as d:
            run = Path(d)
            (run/'phenotype.upload').write_bytes(b'Accession\ttrait\ns1\t1\n')
            script = server.bridge_script(run, 'GPall', 'test', 'test.txt')
            script = script.replace('Option Explicit', 'Option Explicit\nDim crt, fixture, requests\nfixture = String(420123, "X")\nrequests = 0')
            mock = '''Function ExecRemote(command)
  Dim offset, count, parts
  If InStr(command, "wc -c") > 0 Then
    ExecRemote = CStr(Len(fixture))
  Else
    parts = Split(command, "skip=")
    offset = CLng(Split(parts(1), " ")(0))
    parts = Split(command, "count=")
    count = CLng(Split(parts(1), " ")(0))
    If count > 4096 Then Err.Raise 2001
    requests = requests + 1
    ExecRemote = Mid(fixture, offset + 1, count)
  End If
End Function'''
            script = re.sub(r'Function ExecRemote\(command\).*?End Function', lambda m: mock, script, flags=re.S)
            script += '''
Call SyncLog(False)
If logOffset <> 16384 Or requests <> 4 Then Err.Raise 2002
Call SyncLog(True)
If logText <> fixture Or logOffset <> Len(fixture) Then Err.Raise 2003
Call SyncLog(False)
If requests <> 103 Then Err.Raise 2004
fixture = fixture & "APPENDED"
Call SyncLog(False)
If logText <> fixture Then Err.Raise 2005
'''
            target = run/'log-test.vbs'
            target.write_text('\n'.join(script.splitlines()[2:]), 'utf-16')
            result = subprocess.run(['cscript.exe', '//nologo', str(target)], capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((run/'remote.log').read_text('utf-16'), 'X'*420123+'APPENDED')

    def test_recovery_bridge_compiles(self):
        with tempfile.TemporaryDirectory(dir=server.ROOT) as d:
            run=Path(d)
            (run/'phenotype.upload').write_bytes(b'Accession\ttrait\ns1\t1\n')
            meta=dict(population='GPall', directory='test', filename='test.txt', label='test')
            script=server.recovery_script(run,meta)
            self.assertNotIn('Call Stage("prepare",',script)
            self.assertIn('Call WaitStage("admin")',script)
            self.assertIn('/admin.started',script)
            syntax=run/'recovery-syntax.vbs'
            syntax.write_text('\n'.join(script.replace('Option Explicit','Option Explicit\nDim crt').splitlines()[2:]),'utf-16')
            result=subprocess.run(['cscript.exe','//nologo',str(syntax)],capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_data_source_paths(self):
        for pop in ('GPall', 'GPallInd', 'GPallJap'):
            source = server.data_source({}, pop)
            self.assertIn(pop + '_chr{chr}.tped', source['tped'])
            custom = dict(source, tped='/new/seq_chr{chr}.tped', kinship='/new/kin.kinf', sample='/new/sample.list', maf='/new/chr{chr}.vcf')
            self.assertEqual(server.data_source({'dataSource': custom}, pop), custom)
            for invalid in ('relative_chr{chr}.tped', '/data/../chr{chr}.tped', '/data/chr{chr}.tped;id', '/data/chr01.tped'):
                with self.assertRaises(ValueError):
                    server.data_source({'dataSource': dict(custom, tped=invalid)}, pop)

    def data(self):
        return dict(population='GPall', directory='test_2026', filename='LSH2022_Re-fitting.txt',
                    content=base64.b64encode(b'Accession\ttrait\ns1\t1\n').decode(), logDirectory=str(server.ROOT))

    def test_populations_and_log_name(self):
        for pop in ('GPall', 'GPallInd', 'GPallJap'):
            data = self.data(); data['population'] = pop
            result = server.validate(data)
            self.assertEqual(result[0], pop)
            self.assertEqual(result[4].name, 'LSH2022_Re-fitting_GWASrun.log')

    def test_reject_shell_injection_and_empty_upload(self):
        for field, value in [('directory', '../escape'), ('directory', 'a;rm'), ('population', 'GPall;id'),
                             ('filename', '../input.txt'), ('filename', 'a$(id).txt'), ('content', '')]:
            data = self.data(); data[field] = value
            with self.assertRaises(ValueError): server.validate(data)

    def test_existing_labels_and_duplicate_columns(self):
        self.assertEqual(server.batch_label({}, 'LSH2022_Re-fitting.txt'), 'LSH2022')
        self.assertEqual(server.phenotype_traits(b'Accession\tLSH2022_height\n', 'GPall', 'LSH2022'), ['GPall_LSH2022_height'])
        for header in (b'Accession\tHeight\tHeight\n', b'Accession\tHeight\tLSH2022_Height\n', b'Accession\ta;id\n'):
            with self.assertRaises(ValueError): server.phenotype_traits(header, 'GPall', 'LSH2022')

    def test_vbs_compiles_without_executing_main(self):
        with tempfile.TemporaryDirectory(dir=server.ROOT) as d:
            run = Path(d)
            (run / 'phenotype.upload').write_bytes(b'Accession\ttrait\ns1\t1\n')
            script = server.bridge_script(run, 'GPall', 'test', 'test.txt')
            # CScript has no SecureCRT crt global; declare a dummy. Main is not called.
            script = script.replace('Option Explicit', 'Option Explicit\nDim crt')
            syntax = run / 'syntax.vbs'
            syntax.write_text('\n'.join(script.splitlines()[2:]), 'utf-16')
            result = subprocess.run(['cscript.exe', '//nologo', str(syntax)], capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


class WorkflowTests(unittest.TestCase):
    def case(self, pop, failure='', source_mode=''):
        with tempfile.TemporaryDirectory(prefix='gwas-test-', dir=server.ROOT) as d:
            root = Path(d)
            def write(relative, content):
                p = root / relative; p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(content, 'utf-8', newline='\n')
                return p
            unix = '/' + root.drive[0].lower() + root.as_posix()[2:]
            script = server.WORKFLOW.read_text('utf-8')
            script = script.replace('BASE=/data9/home/yzhao/GWAS_IRGSP1.0', 'BASE=' + server.shlex.quote(unix))
            script = script.replace('/public/home/yzhao/IRGSP1.0_GPall/MAF/', unix + '/MAF/')
            script = script.replace('/data5/home/yzhao/program/compared_peakSNP_with_RiceNaviGene_1.0.pl', unix + '/annotate.pl')
            write('workflow.sh', script)
            source = f'000data_prepare/EMMAx.Data/{pop}_miss20'
            for name in (f'{pop}.sampleList', f'{pop}_chr01.tped', f'{pop}_allChr_snpNonHet.hIBS.kinf'):
                write(source + '/' + name, 'fixture\n')
            for chromosome in range(1,13):
                write(f'MAF/{pop}_chr{chromosome:02d}_MAF.vcf', 'fixture\n')
                write(source+f'/{pop}_chr{chromosome:02d}.tped', 'fixture\n')
            program = '000data_prepare/program/'
            write(program + 'select.pheno.only.pl', 'open(my $in, "<", $ARGV[1]) or die $!; open(my $f, ">", $ARGV[2]) or die $!; while(<$in>) {s/\\t/,/g; print $f $_;}')
            maximum = 11 if failure == 'chromosome' else 12
            converter = 'my $p=$ARGV[1]; my @n=("${p}_height.pheno", "${p}_height.tfam"); for my $i (1..MAXIMUM) {push @n,"${p}_height_chr${i}.tped"; push @n,"${p}_height_chr${i}.tfam";} for my $n (@n) {open(my $f, ">", "EMMAx.Data/$n") or die $!; print $f "fixture\\n";}'
            write(program + 'EMMAx-Step2.ChangetoEMMAxPhenotype.tfam.pl', converter.replace('MAXIMUM', str(maximum)))
            if not failure:
                target = root / (program+'EMMAx-Step2.ChangetoEMMAxPhenotype.tfam.pl')
                target.write_bytes((server.ROOT/'provided_program'/target.name).read_bytes().replace(b'\r\n', b'\n'))
            write(program + 'qsub_nodeAdmin.pl', 'my ($n,$c)=@ARGV; open(my $script, ">", $n) or die $!; print $script "$c\\n"; close $script; my $rc=system("$c >$n.o123 2>$n.e123"); print "123.mock\\n";')
            worker = 'open(my $f, ">", "EMMAx.Result/hIBS/$ARGV[0].hIBS") or die $!; print $f "CHROM POS P\\n01 100 6\\n01 2000000 2\\n"; for my $ext ("log","ps","reml") {open(my $g, ">", "EMMAx.Result/hIBS/output/$ARGV[0]_hIBS.$ext") or die $!; print $g "fixture\\n";}'
            if failure == 'worker': worker = 'exit 17;'
            if failure == 'errorlog': worker += 'print STDERR "ERROR intentional fixture\\n";'
            if failure == 'missing': worker = 'exit 0;'
            write(program + 'EMMAx.run-get.pValue.hIBS.pl', worker)
            write(program + 'plot_Allpicture_ofOneTrait.pl', 'open(my $f, ">", "plot/$ARGV[1].png") or die $!; print $f "fixture\\n";')
            write(program + 'plot_manhattan.pl', 'open(my $f, ">", "chromosome/$ARGV[1].png") or die $!; print $f "fixture\\n";')
            write(program + 'plot_manhattan_eachChr.pl', 'for my $c (1..12) {my $n=sprintf("%02d",$c); open(my $f, ">", "plot/$ARGV[1]_chr$n.png") or die $!; print $f "fixture\\n";}')
            write(program + 'plot_QQplot2023.pl', 'for my $n ("$ARGV[1]_QQplot.png", "$ARGV[1].value") {open(my $f, ">", "QQplot/$n") or die $!; print $f "fixture\\n";}')
            if failure == 'plot': write(program+'plot_manhattan.pl', 'exit 23;')
            if failure == 'plotmissing': write(program+'plot_QQplot2023.pl', 'exit 0;')
            write(program + 'stat.allSNPSignificantPos.inOneDir.list.pl', 'open(my $f, ">", "sigSNP/$ARGV[2].list") or die $!; for my $n (glob("$ARGV[0]/*.hIBS")) {print $f "$n\\n";}')
            if not failure:
                target = root/(program+'stat.allSNPSignificantPos.inOneDir.list.pl')
                target.write_bytes((server.ROOT/'provided_program'/target.name).read_bytes().replace(b'\r\n', b'\n'))
            write('annotate.pl', 'open(my $f, ">", $ARGV[1]) or die $!; print $f "gene\\n";')
            write('bin/id', '#!/bin/sh\necho yzhao\n')
            queue_counter = server.shlex.quote(unix+'/queue-counter')
            queue_script = '#!/bin/sh\nn=$(cat '+queue_counter+' 2>/dev/null || echo 0)\nn=$((n+1))\necho "$n" > '+queue_counter+'\ns=C\nif [ "$n" = 1 ]; then s=R; fi\necho "123.mock fixture yzhao 00:00:01 $s low"\n'
            write('bin/qstat', '#!/bin/sh\nexit 1\n' if failure == 'qstat' else queue_script)
            write('bin/qsub', '#!/bin/sh\nexit 0\n')
            write('bin/sleep', '#!/bin/sh\necho "$1" >> ' + server.shlex.quote(unix+'/sleep-events') + '\nexit 0\n')
            # Pre-existing manual files and an old error log must be ignored, not moved.
            legacy = ['analysis/EMMAx.Data/'+pop+'_legacy_chr01.tped',
                      'analysis/EMMAx.Result/hIBS/'+pop+'_legacy_chr01.hIBS',
                      'analysis/EMMAx.Result/hIBS/output/'+pop+'_legacy_chr01_hIBS.log',
                      'analysis/EMMAx.Result/hIBS/ManualHeadingDate/old.hIBS',
                      'analysis/old_plotGWAS.sh.o42', 'analysis/sigSNP/'+pop+'.list']
            for name in legacy: write(name, 'ERROR old content must remain unchanged\n')
            env = dict(os.environ)
            # Git Bash maps C:/ paths and provides its own Perl.
            cmd_prefix = 'export PATH=' + server.shlex.quote(unix + '/bin') + ':$PATH; '
            first_batch_files = {}
            for year in (2022, 2023):
                rid = f'fixture{year}'
                stem = f'LSH{year}_Re-fitting'
                label = f'LSH{year}'
                archive = f'{stem}__{rid}'
                if source_mode == 'custom':
                    archive = f'{stem}_20260920-0830'
                    write(f'.gwas-web/{rid}/archive-name.txt', archive)
                write(f'.gwas-web/{rid}/phenotype.upload', f'Accession\t{label}_height\nGP001\t1\n')
                write(f'.gwas-web/{rid}/expected_traits.txt', f'{pop}_{label}_height\n')
                if source_mode == 'custom':
                    for chromosome in range(1,13):
                        write(f'resequence/new_{chromosome:02d}.tped', 'new genotype\n')
                        write(f'resequence/maf_{chromosome:02d}.vcf', 'new maf\n')
                    write('resequence/samples.list', 'new samples\n')
                    write('resequence/kin.kinf', 'new kinship\n')
                    write(f'.gwas-web/{rid}/data-source.txt', '\n'.join(['reseq-v2', unix+'/resequence/samples.list', unix+'/resequence/new_{chr}.tped', unix+'/resequence/kin.kinf', unix+'/resequence/maf_{chr}.vcf'])+'\n')
                if source_mode == 'changed' and year == 2023:
                    write(source+f'/{pop}_chr01.tped', 'changed source bytes\n')
                for stage in ('prepare', 'admin', 'final'):
                    args = ['bash', unix + '/workflow.sh', stage, pop, 'analysis', stem+'.txt', rid, label]
                    command = cmd_prefix + ' '.join(server.shlex.quote(a) for a in args)
                    result = subprocess.run([BASH, '-c', command], env=env, capture_output=True, timeout=90)
                    log = (root / f'.gwas-web/{rid}/run.log').read_text('utf-8', errors='replace')
                    if stage == 'prepare':
                        self.assertEqual((root/f'analysis/phenotype/{stem}.txt').read_bytes(),
                                         (root/f'.gwas-web/{rid}/phenotype.upload').read_bytes())
                        self.assertLess(log.index('UPLOAD_READY:'), log.index('perl program/select.pheno.only.pl'))
                        if source_mode == 'changed' and year == 2023:
                            self.assertEqual((root/f'analysis/tped/{pop}_chr01.tped').read_text(), 'changed source bytes\n')
                        if source_mode == 'custom':
                            self.assertEqual((root/f'analysis/tped/{pop}_chr01.tped').read_text(), 'new genotype\n')
                    if failure and stage == ('prepare' if failure == 'chromosome' else 'admin'):
                        self.assertNotEqual(result.returncode, 0, log)
                        self.assertFalse((root / f'analysis/sigSNP/{archive}').exists())
                        self.assertIn('ERROR:', log)
                        return
                    self.assertEqual(result.returncode, 0, log + str(result.stderr))
                result_dir = root / f'analysis/EMMAx.Result/hIBS/{archive}'
                self.assertEqual(len(list(result_dir.glob('*.hIBS'))), 12)
                self.assertEqual(len(list((result_dir/'output').iterdir())), 36)
                self.assertFalse((root/'analysis/.gwas-active').exists())
                jobs = (root/f'analysis/.gwas-runs/{rid}/emmax.jobs').read_text()
                self.assertEqual(len(jobs.splitlines()), 12)
                self.assertNotIn('legacy', jobs)
                self.assertNotIn(f'LSH{2023 if year==2022 else 2022}', jobs)
                self.assertNotIn('Re-fitting', jobs)  # Existing year tag is preserved once.
                plot_jobs = (root/f'analysis/.gwas-runs/{rid}/plot.jobs').read_text()
                self.assertEqual(len(plot_jobs.splitlines()), 1)
                self.assertEqual((root/'sleep-events').read_text().splitlines().count('5'), (year-2021)*13)
                self.assertIn('120', (root/'sleep-events').read_text().splitlines())
                self.assertEqual(len(list((root/f'analysis/.gwas-runs/{rid}/status').glob('*.queue-C'))), 13)
                if source_mode == 'custom':
                    self.assertEqual((root/f'analysis/phenotype/{pop}.sampleList').read_text(), 'new samples\n')
                peaks = (root/f'analysis/sigSNP/{archive}/{pop}.list').read_text()
                self.assertEqual(len(peaks.splitlines()), 13)  # Real peak program includes header.
                self.assertTrue(all(label in line for line in peaks.splitlines()[1:]))
                self.assertTrue((root/f'analysis/sigSNP/{archive}/{pop}_geneList.csv').exists())
                self.assertEqual(len(list((root/f'analysis/sh.file/emmax/{archive}').glob('*.o*'))), 12)
                self.assertEqual(len(list((root/f'analysis/sh.file/plotGWAS/{archive}').glob('*.o*'))), 1)
                self.assertEqual(len(list((root/f'analysis/sh.file/emmax/{archive}').glob('*.sh'))), 12)
                self.assertEqual(len(list((root/f'analysis/sh.file/plotGWAS/{archive}').glob('*.sh'))), 1)
                self.assertFalse(list((root/'analysis/EMMAx.Data').iterdir()))
                self.assertFalse(list((root/'analysis/tped').iterdir()))
                for name in legacy:
                    if '/EMMAx.Data/' in name:
                        self.assertFalse((root/name).exists())
                        continue
                    self.assertEqual((root/name).read_text(), 'ERROR old content must remain unchanged\n')
                for p, content in first_batch_files.items():
                    self.assertEqual(p.read_bytes(), content, f'Old batch changed: {p}')
                if year == 2022:
                    first_batch_files = {p:p.read_bytes() for p in (root/'analysis').rglob('*')
                                         if p.is_file() and 'LSH2022' in str(p)}

    def test_three_populations(self):
        for pop in ('GPall', 'GPallInd', 'GPallJap'):
            with self.subTest(pop=pop): self.case(pop)

    def test_custom_source_and_changed_source(self):
        self.case('GPall', source_mode='custom')
        self.case('GPall', source_mode='changed')

    def test_stop_on_failures(self):
        for failure in ('worker', 'errorlog', 'missing', 'qstat', 'chromosome', 'plot', 'plotmissing'):
            with self.subTest(failure=failure): self.case('GPallInd', failure)


if __name__ == '__main__': unittest.main(verbosity=2)
