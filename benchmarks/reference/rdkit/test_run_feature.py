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
    def test_smarts_behavior_keeps_ordered_tags_symmetry_and_stereo(self):
        def target(query, id):
            value = reference.smarts_behavior(Chem.MolFromSmarts(query))
            return next(row for row in value['targets'] if row['target'] == id)
        value = target('[#6:2]-[#8:1]', 'smoke:cid_702')
        self.assertEqual(value['matches'], [[1, 2]])
        self.assertEqual(value['tagged_matches'], [[2, 1]])
        self.assertEqual(target('[#6]=[#8]', 'smoke:cid_702')['matches'], [])
        self.assertEqual(len(target('[c:1]1[c:2]cccc1', 'smoke:cid_241')['matches']), 12)
        self.assertEqual(target('C[C@@H](C(=O)O)N', 'smoke:cid_5950')['matches'], [[0,1,2,3,4,5]])
        self.assertEqual(target('C[C@H](C(=O)O)N', 'smoke:cid_5950')['matches'], [])
        value = reference.smarts_behavior(Chem.MolFromSmarts('[#6:0]-[#8:0]'))
        self.assertEqual(value['tags'], [{'tag': 0, 'query_atom': 0}, {'tag': 0, 'query_atom': 1}])

    def test_smarts_reference_enumerates_recursion_completely_and_reports_limit(self):
        # Toy structures test the adapter, never augment the scientific corpus.
        molecule = Chem.MolFromSmiles('.'.join(['C'] * 1002))
        with patch.object(reference, 'smarts_targets', return_value=(
                {'max_matches': 2000}, [('test', molecule)])):
            value = reference.smarts_behavior(Chem.MolFromSmarts('[$([#6])]'))
            self.assertEqual(len(value['targets'][0]['matches']), 1002)
        with patch.object(reference, 'smarts_targets', return_value=(
                {'max_matches': 1000}, [('test', molecule)])):
            with self.assertRaisesRegex(RuntimeError, 'resource limit'):
                reference.smarts_behavior(Chem.MolFromSmarts('*'))

    def test_disconnected_precheck_preserves_rdkit_full_mappings(self):
        params = Chem.SubstructMatchParameters()
        params.useChirality = True
        params.uniquify = False
        params.maxMatches = params.maxRecursiveMatches = 0
        for text in ['O.O', '[$(CO)].O', 'C[C@@H](C(=O)O)N.O']:
            query = Chem.MolFromSmarts(text)
            value = reference.smarts_behavior(query)
            for actual, (_, molecule) in zip(value['targets'], reference.smarts_targets()[1]):
                expected = sorted(map(list, molecule.GetSubstructMatches(query, params)))
                self.assertEqual(actual['matches'], expected)
        # Every target lacks sulfur: a water-fragment factorial search is unnecessary.
        hydrate = Chem.MolFromSmarts('O.O.O.O.O.O.O.O.O.O.[O-]S(=O)(=O)[O-].[Na+].[Na+]')
        self.assertTrue(all(not row['matches'] for row in reference.smarts_behavior(hydrate)['targets']))

    def test_mdl_aromaticity_clears_default_atom_and_bond_flags(self):
        source = 'c1cc[nH]c1.c1ccccc1-c2ccccc2'
        mdl = self.value(source, 'algo.aromaticity.mdl')['records']
        default = self.value(source, 'algo.aromaticity.rdkit-like')['records']
        self.assertEqual(mdl[0]['atom_aromatic'], [False] * 5)
        self.assertTrue(all(not bond['value'] for bond in mdl[0]['bond_aromatic']))
        self.assertEqual(default[0]['atom_aromatic'], [True] * 5)
        self.assertTrue(all(bond['value'] for bond in default[0]['bond_aromatic']))
        self.assertEqual(mdl[1]['atom_aromatic'], [True] * 12)
        self.assertEqual(len(mdl[1]['bond_aromatic']), 13)
        self.assertEqual(sum(not bond['value'] for bond in mdl[1]['bond_aromatic']), 1)
        # The reference operates on a copy and cannot leave its input under MDL.
        molecule = Chem.MolFromSmiles('c1cc[nH]c1')
        reference.aromaticity_record({'record_index': 0, 'title': '', 'mol': molecule}, mdl=True)
        self.assertTrue(all(atom.GetIsAromatic() for atom in molecule.GetAtoms()))
        self.assertTrue(all(bond.GetIsAromatic() for bond in molecule.GetBonds()))

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

    def value(self,text,feature='io.smiles.parse',path='input.smi'):
        result=runner.run({'feature':feature,'inputs':[{'path':path,'text':text}]})['results'][0]
        self.assertEqual(result['status'],'ok',result)
        return result['value']

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


if __name__ == '__main__':
    unittest.main()
