import sys, tempfile, zipfile, json, csv, os, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from qc_tools import qc,trim,bam,report
from prepare_review import fasta
class Tests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.old=os.getcwd();os.chdir(self.tmp.name)
 def tearDown(self):os.chdir(self.old);self.tmp.cleanup()
 def zip(self,name,adapter='pass',quality='pass',n=100):
  with zipfile.ZipFile(name,'w') as z:z.writestr('test/fastqc_data.txt','>>Basic Statistics\tpass\nTotal Sequences\t%d\n>>END_MODULE\n>>Adapter Content\t%s\n>>END_MODULE\n>>Per base sequence quality\t%s\n>>END_MODULE\n'%(n,adapter,quality))
  return name
 def test_qc(self):
  p=[self.zip('r1.zip'),self.zip('r2.zip',quality='warn')]
  q=qc('a','r','trim',p);self.assertTrue(q['adapter_pass']);self.assertTrue(q['quality_pass']);self.assertEqual(q['total_reads'],200)
  self.zip('r2.zip',adapter='fail');self.assertFalse(qc('a','r','raw',p)['adapter_pass'])
  self.zip('r2.zip',n=99);self.assertFalse(qc('a','r','raw',p)['valid'])
 def test_counts(self):
  Path('log').write_text('Input Read Pairs: 100 Both Surviving: 80 Forward Only Surviving: 10 Reverse Only Surviving: 5 Dropped: 5')
  self.assertEqual(trim('a','r','log')['reads_after_trim'],175)
  Path('count').write_text('150');Path('depth').write_text('contig\tpositions\tmean_depth\n1\t100\t10\n2\t300\t2\n')
  self.assertEqual(bam('a','r','run','count','depth')['mean_depth'],4)
 def test_report_waiting(self):
  Path('manifest').write_text('animal_id\trun\na\tr1\na\tr2\nb\tr3\n')
  data=[dict(run='r1',animal='a',kind='run',mapped_primary_reads=100,mean_depth=2),dict(run='b',animal='b',kind='final',mapped_primary_reads=200,mean_depth=4)]
  for i,d in enumerate(data):Path(str(i)+'.json').write_text(json.dumps(d))
  report('manifest',['0.json','1.json'])
  with open('animal_statistics.tsv') as f: rows=list(csv.DictReader(f,delimiter='\t'))
  self.assertEqual(rows[0]['pending_runs'],'r2');self.assertEqual(rows[0]['status'],'WAITING');self.assertEqual(rows[1]['status'],'FINAL')
 def test_adapters(self):
  Path('a.fa').write_text('>custom\nACGTN\n');self.assertEqual(fasta('a.fa'),[('custom','ACGTN')])
  Path('a.fa').write_text('>bad\nNOTDNA\n')
  with self.assertRaises(ValueError):fasta('a.fa')
if __name__=='__main__':unittest.main()
