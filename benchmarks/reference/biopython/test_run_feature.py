from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import run_feature as reference


class BiopythonReferenceTests(unittest.TestCase):
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

    def test_dssp_receives_original_bytes_including_archive_metadata(self):
        text='data_source\n_pdbx_database_related.db_name PDB\n'
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.cif'
            path.write_text(text)
            seen=[]
            def dssp(model, source, **kwargs):
                seen.append((Path(source),Path(source).read_text()))
                return {}
            parser=lambda **kwargs:SimpleNamespace(get_structure=lambda *args:[object()])
            with patch.object(reference.subprocess,'run'),patch.object(reference,'dssp_extended_rows',return_value=[]):
                value=reference.dssp_summary(path,parser,lambda *args:{},dssp)
            self.assertEqual(seen,[(path,text)])
            self.assertEqual(value['status'],'no_analyzable_residues')


if __name__=='__main__':
    unittest.main()
