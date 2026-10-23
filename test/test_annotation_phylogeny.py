from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from annotation_phylogeny import load_phylogeny


class PhylogenyTests(unittest.TestCase):
    def test_scientific_names_and_custom_assembly_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tree = root / 'tree.nwk'
            tree.write_text('(Beta species,Alpha species);')
            taxonomy = root / 'taxonomy.tsv'
            taxonomy.write_text('id\tname\nsample-A\tAlpha species\nsample-B\tBeta species\n')
            _, _, order, missing = load_phylogeny(tree, ['sample-A', 'sample-B'], taxonomy)
            self.assertEqual(order, ['sample-B', 'sample-A'])
            self.assertEqual(missing, [])

    def test_multiple_assemblies_of_one_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tree = root / 'tree.nwk'; tree.write_text('(Alpha species);')
            taxonomy = root / 'taxonomy.tsv'
            taxonomy.write_text('id\tspecies\nassembly-1\tAlpha species\nassembly-2\tAlpha species\n')
            parsed, leaves, order, missing = load_phylogeny(tree, ['assembly-2', 'assembly-1'], taxonomy)
            self.assertEqual(order, ['assembly-1', 'assembly-2'])
            self.assertEqual(leaves[id(parsed.get_terminals()[0])], order)
            self.assertEqual(missing, [])

    def test_conventional_ids_without_taxonomy(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree = Path(tmp) / 'tree.nwk'; tree.write_text('(Hanseniaspora uvarum);')
            _, _, order, missing = load_phylogeny(tree, ['hanseniaspora-uvarum'])
            self.assertEqual(order, ['hanseniaspora-uvarum'])
            self.assertEqual(missing, [])


if __name__ == '__main__':
    unittest.main()
