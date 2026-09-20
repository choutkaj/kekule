from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import run_feature as reference


class BiopythonReferenceTests(unittest.TestCase):
    def test_version_probe_is_resource_bounded(self):
        with patch.object(reference.shutil,'which',return_value='mkdssp'),patch.object(reference.subprocess,'run',return_value=SimpleNamespace(stdout='mkdssp version 4.6.1')) as execute:
            self.assertIn('4.6.1',reference.dssp_reference('1.87')['version'])
            self.assertEqual(execute.call_args.kwargs['timeout'],120)

    def evaluate(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.cif'
            path.write_text(text)
            return reference.evaluate('io.mmcif.parse',path,reference.import_biopython())

    def test_all_categories_and_missing_tokens_are_preserved(self):
        value=self.evaluate('data_example\n_custom.first ?\n_custom.second .\nloop_\n_atom_site.id\n_atom_site.label_entity_id\n_atom_site.pdbx_formal_charge\n1 9 2\n')
        self.assertEqual(value['blocks'][0]['values'],{
            '_custom.first':['?'],'_custom.second':['.'],'_atom_site.id':['1'],
            '_atom_site.label_entity_id':['9'],'_atom_site.pdbx_formal_charge':['2']})

    def test_multiple_blocks_cannot_silently_merge(self):
        with self.assertRaises((ValueError,KeyError)):
            self.evaluate('data_first\n_x.a 1\ndata_second\n_x.a 2\n')

    def test_embedded_quotes_follow_cif_11_delimiter_rules(self):
        value = self.evaluate('data_x\n_x.text \'don\'t stop\'\nloop_\n_x.name\n_x.id\n"a"b" 1\n')
        self.assertEqual(value['blocks'][0]['values'], {
            '_x.text': ["don't stop"], '_x.name': ['a"b'], '_x.id': ['1']})
        for literal in ["'value'junk", '"value"junk', "'value'#comment"]:
            with self.subTest(literal=literal), self.assertRaises(ValueError):
                self.evaluate(f'data_x\n_x.text {literal}\n')

    def test_multiline_trailing_whitespace_is_biopython_policy(self):
        value = self.evaluate('data_x\n_x.text\n;  first  \n  second\t\n;\n_x.quoted \' padded \'\n')
        self.assertEqual(value['blocks'][0]['values'], {
            '_x.text': ['  first\n  second'], '_x.quoted': [' padded ']})

    def test_dssp_receives_original_bytes_including_archive_metadata(self):
        text='data_source\n_pdbx_database_related.db_name PDB\n'
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.cif'
            path.write_text(text)
            seen=[]
            def execute(command, **kwargs):
                source=Path(command[3])
                seen.append((source,source.read_text()))
                self.assertEqual(kwargs['timeout'],120)
                return SimpleNamespace(stdout='')
            parser=lambda **kwargs:SimpleNamespace(get_structure=lambda *args:[object()])
            with patch.object(reference.subprocess,'run',side_effect=execute),patch.object(reference,'dssp_extended_rows',return_value=[]):
                value=reference.dssp_summary(path,parser,lambda *args:{},lambda *args,**kwargs:{})
            self.assertEqual(seen,[(path,text),(path,text)])
            self.assertEqual(value['status'],'no_analyzable_residues')

    def test_beta_partner_identity_comes_from_explicit_columns(self):
        row=list(' '*80)
        row[:5]='    1'
        row[9]='1'
        row[25:29]='   7'
        row[29:33]='  13'
        text='  #  RESIDUE\n'+''.join(row)
        self.assertEqual(reference.beta_partner_indices(text),{1:[7,13]})
        with self.assertRaises(ValueError):
            reference.beta_partner_indices(text+'\n'+''.join(row))
        self.assertEqual(reference.partner_identity(('A',(' ',13,'B'))),
                         {'partner_chain_id':'A','partner_sequence_id':13,'partner_insertion_code':'B'})

    def test_omega_uses_four_backbone_atoms_and_stops_at_breaks(self):
        from Bio.PDB.vectors import Vector
        class Residue(dict):
            def __init__(self,index,atoms):
                super().__init__({name:SimpleNamespace(get_vector=lambda xyz=xyz:Vector(*xyz)) for name,xyz in atoms.items()})
                self.id=(' ',index,' ')
        keys=[('A',(' ',1,' ')),('A',(' ',2,' '))]
        model={'A':[Residue(1,{'CA':(0,1,0),'C':(0,0,0)}),Residue(2,{'N':(1,0,0),'CA':(1,0,1)})]}
        assignments={keys[0]:[1],keys[1]:[2]}
        self.assertAlmostEqual(abs(reference.omega_angle(model,keys,assignments,0)),90.0)
        self.assertIsNone(reference.omega_angle(model,keys,assignments,1))
        assignments[keys[1]]=[3]
        self.assertIsNone(reference.omega_angle(model,keys,assignments,0))


if __name__=='__main__':
    unittest.main()
