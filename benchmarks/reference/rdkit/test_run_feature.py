import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
from rdkit import Chem
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from reference.rdkit import run_feature as reference
from reference.rdkit import molecule as strict

spec = importlib.util.spec_from_file_location('benchmark_runner', Path(__file__).parents[1] / 'run.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class StrictReferenceTests(unittest.TestCase):
    def test_valence_prepares_oxohalogens_before_property_cache(self):
        for halogen in ['Cl', 'Br', 'I']:
            for oxo_count in [1, 2, 3]:
                source = 'O' + halogen + '(=O)' * oxo_count
                with self.subTest(source=source):
                    atoms = self.value(source, 'algo.valence.rdkit-like')['records'][0]['atoms']
                    self.assertEqual(atoms[1]['formal_charge'], oxo_count)
                    self.assertEqual(atoms[1]['explicit_valence'], oxo_count + 1)
                    self.assertEqual(atoms[0]['implicit_hydrogens'], 1)
                    for oxygen in atoms[2:]:
                        self.assertEqual(oxygen['formal_charge'], -1)
                        self.assertEqual(oxygen['explicit_valence'], 1)

    def test_valence_cleanup_is_independent_and_does_not_sanitize(self):
        # Exercise RDKit's cleanup beyond halogens, without inferring radicals
        # or running aromaticity/stereo perception as part of this feature.
        for source, center in [('CN(=O)=O', 1), ('CC=P(=O)C', 2)]:
            with self.subTest(source=source):
                atoms = self.value(source, 'algo.valence.rdkit-like')['records'][0]['atoms']
                self.assertEqual(atoms[center]['formal_charge'], 1)
                self.assertEqual(atoms[center]['explicit_valence'], 4)
        original = Chem.MolFromSmiles('[CH3]', sanitize=False)
        with patch.object(Chem, 'SanitizeMol', side_effect=AssertionError('full sanitization')):
            result = reference.valence_record({'mol': original, 'record_index': 0, 'title': ''})
        self.assertEqual(result['status'], 'ok')
        self.assertEqual(original.GetAtomWithIdx(0).GetNumRadicalElectrons(), 0)

    def test_valence_cleanup_preserves_source_and_respects_halogen_scope(self):
        for source in ['OCl(=O)=O', 'COCl=O', 'CCl(=O)=O', 'O[Cl+](=O)O']:
            original = Chem.MolFromSmiles(source, sanitize=False)
            original.UpdatePropertyCache(strict=False)
            before = Chem.MolToMolBlock(original, kekulize=False)
            result = reference.valence_record({'mol': original, 'record_index': 0, 'title': ''})
            with self.subTest(source=source):
                self.assertEqual(Chem.MolToMolBlock(original, kekulize=False), before)
                if source == 'COCl=O':
                    self.assertEqual(result['atoms'][2]['formal_charge'], 1)
                elif source in ['CCl(=O)=O', 'O[Cl+](=O)O']:
                    self.assertEqual(result['atoms'],
                                     [reference.valence_atom_json(atom) for atom in original.GetAtoms()])

    def test_sdf_feature_reads_all_records_independently_of_suffix(self):
        first = Chem.MolFromSmiles('C')
        first.SetProp('_Name', 'first')
        second = Chem.MolFromSmiles('O')
        second.SetProp('_Name', 'second')
        text = (Chem.MolToMolBlock(first) + '>  <_Private>\nkept\n\n$$$$\n'
                + Chem.MolToMolBlock(second) + '$$$$\n')
        expected = self.value(text, 'io.sdf.parse', 'input.sdf')
        self.assertEqual([r['title'] for r in expected['records']], ['first', 'second'])
        self.assertEqual(expected['records'][0]['properties'],
                         [{'name': '_Private', 'value': 'kept'}])
        for path in ['input.mol', 'input.mdl']:
            with self.subTest(path=path):
                self.assertEqual(self.value(text, 'io.sdf.parse', path), expected)

    def test_aromatic_observations_preserve_triple_bond_order(self):
        for source in ['C1=CC=CC#C1', 'C1=C(C(=CC#C1)Cl)Cl']:
            with self.subTest(source=source):
                bonds = self.value(source)['records'][0]['components'][0]['bonds']
                triples = [bond for bond in bonds if bond['bond_type'] == 'TRIPLE']
                self.assertEqual(len(triples), 1)
                self.assertTrue(triples[0]['is_aromatic'])
                self.assertEqual(sum(bond['bond_type'] == 'AROMATIC' for bond in bonds), 5)

    def test_writer_version_contract_and_extensions(self):
        mol = Chem.MolFromSmiles('CCO')
        mol.AddConformer(Chem.Conformer(3))
        v2, v3 = Chem.MolToMolBlock(mol), Chem.MolToV3KMolBlock(mol)
        for feature, text, path in [
            ('io.mol.v2000.write', v3, 'output.mol'),
            ('io.mol.v3000.write', v2, 'output.mol'),
            ('io.sdf.v2000.write', v3 + '$$$$\n', 'output.sdf'),
        ]:
            with self.subTest(feature=feature), self.assertRaisesRegex(ValueError, 'must emit'):
                strict.read_written(feature, {'written':[{'path':path,'text':text}]}, runner.MemoryInput)
        for text in ['F[C@H](Cl)Br |&1:1|', 'F[C@H](Cl)Br |&1:1| sample']:
            value = strict.writer_value('io.smiles.isomeric', self.source(text))['records'][0]
            self.assertIn('&', value['identity'])
            self.assertEqual(value['title'], 'sample' if text.endswith('sample') else '')

    def source(self, text, path='input.smi'):
        return runner.MemoryInput(path,text)

    def test_mol_writer_reader_rejects_trailing_records_and_sdf_data(self):
        mol = Chem.MolFromSmiles('CCO')
        for feature, block in [
            ('io.mol.v2000.write', Chem.MolToMolBlock(mol)),
            ('io.mol.v3000.write', Chem.MolToV3KMolBlock(mol)),
        ]:
            def read(text):
                return strict.read_written(feature,
                    {'written': [{'path': 'output.mol', 'text': text}]},
                    runner.MemoryInput)
            with self.subTest(feature=feature):
                self.assertEqual(read(block), read(block + ' \n\t\n'))
                for trailing in ['$$$$\n', '>  <EXTRA>\nvalue\n\n', block]:
                    with self.subTest(trailing=trailing[:20]):
                        with self.assertRaisesRegex(ValueError, 'exactly one complete MOL record'):
                            read(block + trailing)

    def value(self,text,feature='io.smiles.parse',path='input.smi'):
        result=runner.run({'feature':feature,'inputs':[{'path':path,'text':text}]})['results'][0]
        self.assertEqual(result['status'],'ok',result)
        return result['value']

    def test_molecular_writer_reader_rejects_query_atoms_and_bonds(self):
        block = Chem.MolToV3KMolBlock(Chem.MolFromSmiles('C'))
        atom_line = next(line for line in block.splitlines() if line.startswith('M  V30 1 C '))
        query = block.replace(atom_line, atom_line + ' HCOUNT=4')
        valence = block.replace(atom_line, atom_line + ' VAL=4')
        self.assertTrue(Chem.MolFromMolBlock(query).GetAtomWithIdx(0).HasQuery())
        ordinary = Chem.MolFromMolBlock(valence)
        self.assertFalse(ordinary.GetAtomWithIdx(0).HasQuery())
        self.assertEqual(ordinary.GetAtomWithIdx(0).GetTotalNumHs(), 4)
        strict.read_written('io.mol.v3000.write',
            {'written': [{'path': 'output.mol', 'text': valence}]}, runner.MemoryInput)
        bond_query = Chem.MolToV3KMolBlock(Chem.MolFromSmiles('CC')).replace(
            'M  V30 1 1 1 2', 'M  V30 1 8 1 2')
        self.assertTrue(Chem.MolFromMolBlock(bond_query).GetBondWithIdx(0).HasQuery())
        for text in [query, bond_query]:
            with self.assertRaisesRegex(ValueError, 'query atoms or bonds'):
                strict.read_written('io.mol.v3000.write',
                    {'written': [{'path': 'output.mol', 'text': text}]}, runner.MemoryInput)

    def test_graph_collision_is_rejected(self):
        self.assertNotEqual(self.value('C1CCCCC1'),self.value('C1CC1.C1CC1'))

    def test_stereo_placement_collision_is_rejected(self):
        a=self.value('F[C@H](Cl)C[C@@H](F)Br')
        b=self.value('F[C@@H](Cl)C[C@H](F)Br')
        self.assertNotEqual(a,b)
        self.assertNotEqual(a['records'][0]['components'][0]['stereo'],b['records'][0]['components'][0]['stereo'])

    def test_dative_direction_survives(self):
        a=strict.graph(Chem.MolFromSmiles('[NH3]->[Cu+2]'))
        backward=Chem.MolFromSmiles('[NH3]<-[Cu+2]',sanitize=False)
        backward.UpdatePropertyCache(strict=False)
        b=strict.graph(backward)
        self.assertNotEqual(a['bonds'],b['bonds'])

    def test_molfile_memory_input_bonds_coordinates_and_duplicate_fields(self):
        block=Chem.MolToMolBlock(Chem.MolFromSmiles('CCO'))
        a=self.value(block,'io.mol.parse','input.mol')
        bond=block.replace('  1  2  1','  1  2  2',1)
        coordinate=block.replace('    0.0000','   10.0000',1)
        self.assertNotEqual(block,bond)
        self.assertNotEqual(block,coordinate)
        self.assertNotEqual(a,self.value(bond,'io.mol.parse','input.mol'))
        self.assertNotEqual(a,self.value(coordinate,'io.mol.parse','input.mol'))
        sdf=block+'>  <ID>\nfirst\n\n>  <SECOND>\nlast\n\n$$$$\n'
        fields=self.value(sdf,'io.sdf.parse','input.sdf')['records'][0]['properties']
        self.assertEqual(fields,[{'name':'ID','value':'first'},{'name':'SECOND','value':'last'}])
        self.assertNotEqual(self.value(sdf,'io.sdf.parse','input.sdf'),self.value(sdf.replace('first','changed'),'io.sdf.parse','input.sdf'))
        duplicate=sdf.replace('<SECOND>', '<ID>')
        result=runner.run({'feature':'io.sdf.parse','inputs':[{'path':'input.sdf','text':duplicate}]})['results'][0]
        self.assertEqual(result['status'],'error')
        self.assertIn('duplicate SDF', result['message'])

    def test_sdf_private_source_fields_survive_without_internal_properties(self):
        block = Chem.MolToMolBlock(Chem.MolFromSmiles('CCO'))
        sdf = block + '>  <_Private>\nfirst\nsecond\n\n>  <ID>\nlast\n\n$$$$\n'
        expected = [{'name': '_Private', 'value': 'first\nsecond'},
                    {'name': 'ID', 'value': 'last'}]
        for feature in ['io.mol.parse', 'io.sdf.parse', 'stereo.representation', 'stereo.perception']:
            with self.subTest(feature=feature):
                fields = self.value(sdf, feature, 'input.sdf')['records'][0]['properties']
                self.assertEqual(fields, expected)
        duplicate = sdf.replace('<ID>', '<_Private>')
        with self.assertRaisesRegex(ValueError, 'duplicate SDF'):
            strict.sdf_record(duplicate)

    def test_algorithm_adapters_use_native_results(self):
        from rdkit.Chem import Descriptors, rdMolDescriptors
        value=self.value('C1CCCCC1','algo.rings.fast')['records'][0]
        self.assertEqual(value['bond_in_ring'][0],{'begin_atom_index':0,'end_atom_index':1,'value':True})
        for text in ['CCCC', '[H]C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]']:
            value=self.value(text,'descriptor.rotatable-bonds.rdkit-strict')['records'][0]
            params=Chem.SmilesParserParams(); params.removeHs=False
            mol=Chem.MolFromSmiles(text,params)
            self.assertEqual(value['count'], rdMolDescriptors.CalcNumRotatableBonds(mol,rdMolDescriptors.NumRotatableBondsOptions.Strict))
        value=self.value('[NH4+]','descriptor.molecular')['records'][0]
        self.assertEqual(value['average_mass_da'],Descriptors.MolWt(Chem.MolFromSmiles('[NH4+]')))
        value=self.value('[#6]-[#8]  external query','query.smarts')['records'][0]
        self.assertEqual((value['smarts'],value['title'],value['atom_count']),('[#6]-[#8]','external query',2))

    def test_hydrogen_removal_uses_unmodified_rdkit_defaults(self):
        text='[H:7]C'
        params=Chem.SmilesParserParams(); params.removeHs=False
        mol=Chem.MolFromSmiles(text,params)
        expected=strict.graph(Chem.RemoveHs(Chem.AddHs(mol)))
        value=self.value(text,'chem.hydrogen-transforms')['records'][0]['round_trip']
        self.assertEqual(value,expected)

    def test_fast_ring_membership_does_not_require_selected_ring_perception(self):
        # A bridge joins two rings; its atoms can be cyclic while its bond is not.
        mol = Chem.MolFromSmiles('C1CC1C2CC2', sanitize=False)
        record = {'record_index': 0, 'title': '', 'status': 'ok', 'mol': mol}
        with patch.object(Chem, 'GetSymmSSSR', side_effect=AssertionError('selected ring perception')):
            value = reference.ring_record(record)
        self.assertEqual(value['atom_in_ring'], [True] * 6)
        flags = {(b['begin_atom_index'], b['end_atom_index']): b['value']
                 for b in value['bond_in_ring']}
        self.assertFalse(flags[(2, 3)])
        self.assertEqual(sum(flags.values()), 6)
        for order in [Chem.BondType.ZERO, Chem.BondType.DATIVE]:
            mol = Chem.MolFromSmiles('C1CC1', sanitize=False)
            mol.GetBondBetweenAtoms(0, 2).SetBondType(order)
            record['mol'] = mol
            with self.subTest(order=order):
                value = reference.ring_record(record)
                self.assertEqual(value['atom_in_ring'], [False] * 3)
                self.assertFalse(any(b['value'] for b in value['bond_in_ring']))

    def test_writer_identity_rejects_wrong_connectivity_and_stereo(self):
        for a,b in [('C1CCCCC1','C1CC1.C1CC1'),('F[C@H](Cl)C[C@@H](F)Br','F[C@@H](Cl)C[C@H](F)Br')]:
            expected=strict.writer_value('io.smiles.isomeric',self.source(a))
            actual=strict.read_written('io.smiles.isomeric',{'written':[{'path':'output.smi','text':b}]},runner.MemoryInput)
            self.assertNotEqual(expected,actual)

    def test_all_input_categories_and_errors_are_retained(self):
        texts=['F[C@H](Cl)Br','F/C=C/F','*','C(','[H]','CC']
        for feature in ['io.smiles.parse','io.smiles.write','stereo.cip','stereo.representation','stereo.perception']:
            results=runner.run({'feature':feature,'inputs':[{'path':'input.smi','text':text} for text in texts]})['results']
            self.assertEqual(len(results),len(texts))
            self.assertEqual(results[3]['status'],'error')
            self.assertEqual(results[-1]['status'],'ok',results)

    def test_disconnected_algorithms_evaluate_all_components(self):
        for feature in ['algo.rings.fast','descriptor.molecular','algo.substructure.vf2']:
            result=self.value('CC.O',feature)
            self.assertEqual(len(result['records']),2)

    def test_query_mapping_and_unbounded_match_count(self):
        # More than the old 1000-match cap. Toy structures are confined to tests.
        molecule=Chem.MolFromSmiles('C'*1002)
        record={'record_index':0,'title':'','status':'ok','mol':molecule}
        value=reference.substructure_record(record,Chem)
        self.assertEqual(len(value['queries'][0]['matches']),1002)
        carbon_bond=next(q for q in value['queries'] if q['smarts']=='C!@C')
        self.assertIn([0,1],carbon_bond['matches'])
        self.assertIn([1,0],carbon_bond['matches'])

    def test_tetrahedral_parity_has_a_fixed_carrier_convention(self):
        value=self.value('F[C@H](Cl)Br')['records'][0]['components'][0]['stereo']
        self.assertEqual(value,[{'type':'tetrahedral','focus':[1],'carriers':[-1,0,2,3],'parity':1}])

    def test_tetrahedral_convention_follows_positive_signed_volume(self):
        molecule=Chem.MolFromSmiles('C(F)(Cl)(Br)I')
        conformer=Chem.Conformer(5)
        for index,point in enumerate([(0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)]):
            conformer.SetAtomPosition(index,point)
        conformer.Set3D(True)
        molecule.AddConformer(conformer)
        Chem.AssignAtomChiralTagsFromStructure(molecule)
        # det(p0-p3,p1-p3,p2-p3)>0, the public Kekule Clockwise convention.
        self.assertEqual(strict.stereo(molecule)[0][0]['parity'],0)

    def test_v3000_inputs_are_read_directly(self):
        block=Chem.MolToV3KMolBlock(Chem.MolFromSmiles('CCO'))
        value=self.value(block,'io.mol.parse','input.mol')
        self.assertEqual(value['records'][0]['components'][0]['atom_count'],3)

    def test_cip_keeps_explicit_hydrogen_vertices_and_empty_descriptors(self):
        value=self.value('[H]C([H])([H])[H]','stereo.cip')
        self.assertEqual(value['records'][0]['atom_count'],5)
        self.assertEqual(value['records'][0]['atom_descriptors'],[])

    def test_cip_uses_component_order_for_interleaved_source_atoms(self):
        source=Chem.MolFromSmiles('F[C@H](Br)I.Cl')
        interleaved=Chem.RenumberAtoms(source,[0,4,1,2,3])
        record={'record_index':0,'title':'','status':'ok','mol':interleaved}
        value=reference.stereo_cip_record(record,Chem)
        self.assertEqual(value['atom_count'],5)
        self.assertEqual(value['atom_descriptors'][0]['atom_index'],1)


class SdfFramingTests(unittest.TestCase):
    def read(self, text):
        return strict.read_written('io.sdf.v2000.write',
            {'written': [{'path': 'output.sdf', 'text': text}]},
            runner.MemoryInput)['records']

    def test_sdf_literal_delimiters_are_metadata(self):
        for title, name, value in [
            ('price $$$$ sample', 'NOTE', 'normal'),
            ('', 'PRICE$$$$', 'normal'),
            ('', 'NOTE', 'price $$$$ sample'),
        ]:
            with self.subTest(title=title, name=name, value=value):
                mol = Chem.MolFromSmiles('C')
                mol.SetProp('_Name', title)
                block = Chem.MolToMolBlock(mol) + f'>  <{name}>\n{value}\n\n$$$$\n'
                records = self.read(block + Chem.MolToMolBlock(Chem.MolFromSmiles('O')) + '$$$$\n')
                self.assertEqual(len(records), 2)
                self.assertEqual(records[0]['title'], title)
                self.assertEqual(records[0]['properties'], [{'name': name, 'value': value}])
                self.assertEqual(records[1]['title'], '')
                self.assertEqual(records[1]['properties'], [])

    def test_sdf_field_headers_are_only_read_outside_values(self):
        mol = Chem.MolFromSmiles('C')
        mol.SetProp('_Name', '>  <TITLE>')
        value = '>  <NOTE>\nlast line'
        text = (Chem.MolToMolBlock(mol) + f'>  <NOTE>\n{value}\n\n'
                + '>  <_Private>\nkept\n\n$$$$\n')
        record = self.read(text)[0]
        self.assertEqual(record['title'], '>  <TITLE>')
        self.assertEqual(record['properties'], [
            {'name': 'NOTE', 'value': value}, {'name': '_Private', 'value': 'kept'}])

    def test_sdf_writer_reader_requires_complete_record_framing(self):
        block = Chem.MolToMolBlock(Chem.MolFromSmiles('C'))
        for text in [
            block + 'junk\n$$$$\n',
            block + '>  <NOTE>\nvalue $$$$',
            block + '$$$$\n' + block,
            block + '$$$$\n$$$$\n',
            block.replace('M  END\n', '') + '$$$$\n',
        ]:
            with self.subTest(text=text[-60:]), self.assertRaises(ValueError):
                self.read(text)
        self.assertEqual(self.read(block + '$$$$\n'), self.read(block + '$$$$\n \n\t\n'))


if __name__ == '__main__':
    unittest.main()
