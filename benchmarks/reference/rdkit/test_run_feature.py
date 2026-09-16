import importlib.util
from pathlib import Path
import unittest
from rdkit import Chem
import run_feature as reference
import strict

spec = importlib.util.spec_from_file_location('benchmark_runner', Path(__file__).parents[1] / 'run.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class StrictReferenceTests(unittest.TestCase):
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
        a=strict.graph(Chem.MolFromSmiles('[NH3]->[Cu+2]'),reference)
        backward=Chem.MolFromSmiles('[NH3]<-[Cu+2]',sanitize=False)
        backward.UpdatePropertyCache(strict=False)
        b=strict.graph(backward,reference)
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
        expected=strict.graph(Chem.RemoveHs(Chem.AddHs(mol)),reference)
        value=self.value(text,'chem.hydrogen-transforms')['records'][0]['round_trip']
        self.assertEqual(value,expected)

    def test_writer_identity_rejects_wrong_connectivity_and_stereo(self):
        for a,b in [('C1CCCCC1','C1CC1.C1CC1'),('F[C@H](Cl)C[C@@H](F)Br','F[C@@H](Cl)C[C@H](F)Br')]:
            expected=strict.writer_value('io.smiles.isomeric',self.source(a),reference)
            actual=strict.read_written('io.smiles.isomeric',{'written':[{'path':'output.smi','text':b}]},reference,runner.MemoryInput)
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
